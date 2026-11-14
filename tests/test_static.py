"""
Static checks: shell and Python syntax, YAML validity, and the two Nextflow source invariants
behind historical bugs (panel-of-normals handling, output paths under projectDir).
"""
import glob
import os
import py_compile
import re
import subprocess

import pytest

from conftest import REPO


def _files(*patterns):
    out = []
    for p in patterns:
        out += glob.glob(os.path.join(REPO, p), recursive=True)
    return sorted(f for f in out if "/work/" not in f and "/.nextflow" not in f and "/_deps/" not in f)


@pytest.mark.parametrize("path", _files("*.sh", "bin/*.sh"))
def test_shell_syntax(path):
    p = subprocess.run(["bash", "-n", path], capture_output=True, text=True)
    assert p.returncode == 0, p.stderr


@pytest.mark.parametrize("path", _files("*.py", "bin/*.py", "tests/*.py"))
def test_python_syntax(path):
    py_compile.compile(path, doraise=True)


@pytest.mark.parametrize("path", _files("yml/*.yml", ".github/workflows/*.yml", "tests/*.yml"))
def test_yaml_parses(path):
    yaml = pytest.importorskip("yaml")
    with open(path) as fh:
        doc = yaml.safe_load(fh)
    assert isinstance(doc, dict), f"{path}: top level is not a mapping"


def test_mutect2_module_adds_pon_only_through_the_conditional():
    src = open(os.path.join(REPO, "Modules", "Mutect2", "MUTECT2_CALLING.nf")).read()
    assert re.search(r'def pon_arg = params\.mutect2_pon \? "-pon \$\{params\.mutect2_pon\}" : ""', src)
    assert src.count("-pon") == 1, "-pon must appear only in the pon_arg definition"
    assert "${pon_arg}" in src


def test_main_nf_pon_is_optional_and_rn7_has_none():
    src = open(os.path.join(REPO, "main.nf")).read()
    assert re.search(r"params\.mutect2_pon = \(params\.genomes\[params\.genome\]\.mutect2_pon \?", src)
    rn7 = src[src.index("'RN7': ["):]
    assert re.search(r'mutect2_pon:\s*""', rn7[:rn7.index("]")])


def test_no_output_path_is_anchored_to_projectdir():
    src = open(os.path.join(REPO, "main.nf")).read()
    bad = [l for l in src.splitlines() if re.search(r'params\.\w+\s*=\s*"\$projectDir/(RESULTS|pipeline\.log|mkdup_tmp)', l)]
    assert not bad, bad
    assert 'params.outdir = "${launchDir}/RESULTS"' in src


def test_tool_selection_uses_exact_membership_only():
    src = open(os.path.join(REPO, "main.nf")).read()
    assert "params.tool.contains(" not in src and "params.tool.toString().contains(" not in src
    assert "def SELECTABLE_TOOLS = ['ascat', 'manta', 'delly', 'cnvkit']" in src
    assert "manta_bed" not in src                      # stale GRCh38-only Manta parameter (audit item 14)


def test_every_java_step_bounds_its_heap_and_gc_threads():
    # GATK, Picard, SAGE, Conpair (GATK3) and FastQC run on the JVM. Each must cap the heap from task.memory
    # and the garbage-collector threads from task.cpus; otherwise the JVM sizes both from the whole node.
    java_steps = []
    for path in _files("Modules/*.nf", "Modules/*/*.nf"):
        src = open(path).read()
        if re.search(r"^\s*(def cmd = \")?(\$\{params\.gatk\}|picard |java |fastqc )", src, re.M):
            java_steps.append(path)
            assert "ParallelGCThreads=${task.cpus}" in src, f"{os.path.relpath(path, REPO)}: no GC thread bound"
            assert "ConcGCThreads=" in src, f"{os.path.relpath(path, REPO)}: no ConcGCThreads bound"
            if "fastqc " not in src:                         # FastQC sizes its own heap per thread
                assert "-Xmx${xmx}M" in src, f"{os.path.relpath(path, REPO)}: no -Xmx from task.memory"
    assert len(java_steps) >= 15, java_steps


def test_recalibrated_bam_is_merged_once_and_not_resorted():
    # The 20 whole-contig ApplyBQSR chunks are coordinate-sorted, so samtools merge yields the final
    # BAM; a second full sort was removed. -c/-p keep one @RG/@PG instead of per-chunk renamed copies.
    assert not os.path.exists(os.path.join(REPO, "Modules", "RECALIBRATE_SortBam.nf"))
    merge = open(os.path.join(REPO, "Modules", "RECALIBRATE_MergeBam.nf")).read()
    assert "samtools merge -@ ${task.cpus} -c -p -o" in merge and "samtools sort" not in merge
    assert "samtools idxstats" in merge                      # the read-count check moved here
    main = open(os.path.join(REPO, "main.nf")).read()
    assert "RECALIBRATE_SortBam" not in main
    assert main.count("RECALIBRATE_MergeBam(MergeBam_input).set{ RECALIBRATE_out }") == 2


def test_markduplicates_spills_to_the_task_directory_by_default():
    mkdup = open(os.path.join(REPO, "Modules", "MKDUP.nf")).read()
    assert "def tmp = params.mkdup_temp_dir ?: '$PWD/tmp'" in mkdup
    assert "TMP_DIR=${tmp}" in mkdup and "-Djava.io.tmpdir=${tmp}" in mkdup
    main = open(os.path.join(REPO, "main.nf")).read()
    assert "params.mkdup_temp_dir = null" in main and "mkdup_tmp" not in main


def test_consensus_script_defaults_to_primary_chromosomes_only():
    src = open(os.path.join(REPO, "consensus_scenarios.sh")).read()
    assert 'printf "${prefix}%s\\n" $(seq 1 22) X Y' in src
