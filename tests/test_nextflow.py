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


def test_outputs_go_to_the_launch_directory(nf):
    launch = nf["root"] / "elsewhere"
    p, out = preview(nf, "--type", "genome", "--genome", "GRCh38", "--first_step", "mapping", launch=launch)
    assert p.returncode == 0, out
    assert f"Output Dir     : {launch}/RESULTS" in out
    assert (launch / "pipeline.log").exists()
    assert not os.path.exists(os.path.join(REPO, "pipeline.log")) or os.path.getmtime(os.path.join(REPO, "pipeline.log")) < os.path.getmtime(launch / "pipeline.log")
