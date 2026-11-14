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


@pytest.mark.parametrize("path", _files("*.sh", "bin/*.sh", "Benchmark/**/*.sh"))
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
    assert re.search(r"params\.mutect2_pon = params\.mutect2_pon \?: \(params\.genomes\[params\.genome\]\.mutect2_pon \?", src)
    rn7 = src[src.index("'RN7': ["):]
    assert re.search(r'mutect2_pon:\s*""', rn7[:rn7.index("]")])


def test_no_output_path_is_anchored_to_projectdir():
    src = open(os.path.join(REPO, "main.nf")).read()
    bad = [l for l in src.splitlines() if re.search(r'params\.\w+\s*=\s*"\$projectDir/(RESULTS|pipeline\.log|mkdup_tmp)', l)]
    assert not bad, bad
    assert 'params.outdir = "${launchDir}/RESULTS"' in src


def test_consensus_script_defaults_to_primary_chromosomes_only():
    src = open(os.path.join(REPO, "consensus_scenarios.sh")).read()
    assert 'printf "${prefix}%s\\n" $(seq 1 22) X Y' in src
