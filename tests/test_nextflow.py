"""
Nextflow smoke tests with `-preview`: the script compiles, every module resolves, the startup
validation runs, and the genome-specific panel-of-normals behaviour is what the README states.
No process is executed, so no tool or database is needed beyond a few placeholder files.
"""
import gzip
import os
import shutil
import subprocess

import pytest

from conftest import REPO


@pytest.fixture(scope="module")
def nf(tmp_path_factory):
    if not shutil.which("nextflow"):
        pytest.skip("nextflow not on PATH")
    root = tmp_path_factory.mktemp("nf")
    db = root / "db"
    (db / "Databases" / "GRCh38").mkdir(parents=True)
    (db / "Databases" / "GRCh38" / "MuTect2.PON.5210.vcf.gz").write_bytes(b"")   # existence is checked at startup
    (db / "Databases" / "mm39").mkdir(parents=True)
    (db / "Databases" / "mm39" / "PoN.mm39.vcf.gz").write_bytes(b"")
    for r in ("R1", "R2"):
        with gzip.open(root / f"{r}.fastq.gz", "wt") as fh:
            fh.write("")
    sheet = root / "sample.csv"
    sheet.write_text("patient,sample,status,fastq_1,fastq_2\n"
                     f"P1,A,normal,{root}/R1.fastq.gz,{root}/R2.fastq.gz\n"
                     f"P1,A,tumor,{root}/R1.fastq.gz,{root}/R2.fastq.gz\n")
    bad = root / "bad.csv"
    bad.write_text("patient,sample,status,fastq_1,fastq_2\n"
                   f"P 1;$(x),A,normal,{root}/R1.fastq.gz,{root}/R2.fastq.gz\n"
                   f"P 1;$(x),A,tumor,{root}/R1.fastq.gz,{root}/R2.fastq.gz\n")
    return dict(root=root, db=str(db), sheet=str(sheet), bad=str(bad))


def preview(nf, *args, launch=None):
    launch = launch or nf["root"] / "launch"
    launch.mkdir(exist_ok=True)
    cmd = ["nextflow", "run", os.path.join(REPO, "main.nf"), "-preview", "-profile", "conda",
           "--database_path", nf["db"], "--verify_database", "none", "--sample", nf["sheet"], *args]
    p = subprocess.run(cmd, cwd=launch, capture_output=True, text=True, env=dict(os.environ, NXF_ANSI_LOG="false"))
    out = p.stdout + p.stderr
    log = launch / "pipeline.log"                       # main.nf writes its configuration block here too
    if log.exists():
        out += log.read_text()
    return p, out


def test_grch38_preview_compiles_and_uses_the_panel_of_normals(nf):
    p, out = preview(nf, "--type", "genome", "--genome", "GRCh38", "--first_step", "mapping")
    assert p.returncode == 0, out
    assert "Pipeline       : AlexandrovLab/SomaticBunny" in out
    assert "has no Mutect2 panel of normals" not in out
    assert "SomaticBunny consensus (GRCh38): consensus_scenarios.sh" in out


def test_rn7_is_refused_without_opt_in(nf):
    p, out = preview(nf, "--type", "genome", "--genome", "RN7", "--first_step", "mapping")
    assert p.returncode != 0
    assert "not validated" in out and "--allow_unvalidated_genome" in out


def test_rn7_runs_without_pon_when_opted_in(nf):
    p, out = preview(nf, "--type", "genome", "--genome", "RN7", "--first_step", "mapping", "--allow_unvalidated_genome", "true")
    assert p.returncode == 0, out
    assert "Genome RN7 has no Mutect2 panel of normals: Mutect2 runs without -pon" in out


def test_invalid_first_step_and_bad_sample_ids_are_refused(nf):
    p, out = preview(nf, "--type", "genome", "--genome", "GRCh38", "--first_step", "foo")
    assert p.returncode != 0 and "--first_step must be one of" in out
    p, out = preview(nf, "--type", "genome", "--genome", "GRCh38", "--first_step", "mapping", "--sample", nf["bad"])
    assert p.returncode != 0 and "is not a valid identifier" in out


def test_unknown_and_unsupported_tools_are_refused_at_startup(nf):
    # --tool is validated once, by exact name, before any task: a typo, a partial match and a human-only
    # tool on a mouse genome all stop the run with a message naming the offending tool.
    for bad in ("mutec2", "ascatXYZ", "manta,ascatXYZ"):
        p, out = preview(nf, "--type", "genome", "--genome", "GRCh38", "--first_step", "mapping", "--tool", bad)
        assert p.returncode != 0 and "unknown tool(s) in --tool" in out and bad.split(",")[-1] in out, out
    p, out = preview(nf, "--type", "genome", "--genome", "mm39", "--first_step", "mapping", "--tool", "ascat")
    assert p.returncode != 0 and "not supported for genome mm39: ascat" in out, out


def test_outputs_go_to_the_launch_directory(nf):
    launch = nf["root"] / "elsewhere"
    p, out = preview(nf, "--type", "genome", "--genome", "GRCh38", "--first_step", "mapping", launch=launch)
    assert p.returncode == 0, out
    assert f"Output Dir     : {launch}/RESULTS" in out
    assert (launch / "pipeline.log").exists()
    assert not os.path.exists(os.path.join(REPO, "pipeline.log")) or os.path.getmtime(os.path.join(REPO, "pipeline.log")) < os.path.getmtime(launch / "pipeline.log")


def test_mm39_preview_announces_the_nonhuman_consensus(nf):
    # Audit point 10: a mm39 run must say at startup which consensus it produces.
    p, out = preview(nf, "--type", "genome", "--genome", "mm39", "--first_step", "mapping")
    assert p.returncode == 0, out
    assert "SomaticBunny consensus (mm39): consensus_scenarios_nonhuman.sh" in out
    assert "SNVs need 2 of 3 callers, indels 2 of 2" in out


def test_rn7_preview_announces_the_nonhuman_consensus(nf):
    p, out = preview(nf, "--type", "genome", "--genome", "RN7", "--first_step", "mapping", "--allow_unvalidated_genome", "true")
    assert p.returncode == 0, out
    assert "SomaticBunny consensus (RN7): consensus_scenarios_nonhuman.sh" in out


def test_consensus_is_wired_for_every_genome():
    # The POST block must not be limited to human genomes, and each branch must call its own module.
    src = open(os.path.join(REPO, "main.nf")).read()
    post = src[src.index("// POST - filtering and consensus calling"):]
    head = post[:post.index("{")]
    assert "GRCh38" not in head and "GRCh37" not in head
    assert "POST_SCENARIOS(postInput)" in post and "POST_SCENARIOS_NONHUMAN(postInput)" in post
    nonhuman = open(os.path.join(REPO, "Modules", "POST_scenarios_nonhuman.nf")).read()
    assert "consensus_scenarios_nonhuman.sh" in nonhuman and "sage_vcf" not in nonhuman
