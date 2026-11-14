"""
pytest fixtures shared by the SomaticBunny regression tests.

Layout: tests/ lives in the repository root, next to consensus_scenarios.sh and VAF.py.
Synthetic inputs are generated once per session by make_fixtures.py into a temporary directory;
consensus_scenarios.sh is run on them in a few configurations and the outputs are parsed here.

Environment variables (all optional):
  SB_REPO               repository root (default: parent of this directory)
  SB_DKFZ_DIR           a checkout of DKFZBiasFilter (default: $SB_DB/DKFZBiasFilter, then tests/_deps/DKFZBiasFilter)
  SB_DB                 a SomaticBunny database directory; only DKFZBiasFilter is taken from it
  SB_TOOLS_BIN          colon-separated bin directories put in front of PATH for every pipeline tool the
                        tests launch (bcftools, alleleCounter, ...): the production yml/POST.yml and
                        yml/alleleCounter.yml environments in the CI job "production-envs" (tests/README.md)
  SB_BIASFILTER_PYTHON  interpreter for the DKFZ bias filter (production: the Python 2.7 of yml/POST.yml)
  SB_VAF_PYTHON         interpreter for VAF.py (production: the python3 of yml/alleleCounter.yml)
  SB_PRODUCTION_ENVS=1  also assert that these tools are the versions pinned in yml/ (test_production_envs.py)
"""
import os
import py_compile
import shutil
import subprocess
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import make_fixtures as mf   # noqa: E402

REPO = os.environ.get("SB_REPO") or os.path.dirname(HERE)
CONSENSUS = os.path.join(REPO, "consensus_scenarios.sh")
VAF_PY = os.path.join(REPO, "VAF.py")
TOOLS_BIN = os.environ.get("SB_TOOLS_BIN", "")
VAF_PYTHON = os.environ.get("SB_VAF_PYTHON", sys.executable)


def tool_path():
    return f"{TOOLS_BIN}:{os.environ.get('PATH', '')}" if TOOLS_BIN else os.environ.get("PATH", "")


def tool_env(**extra):
    """Environment for every subprocess that runs a pipeline tool: SB_TOOLS_BIN first on PATH, C locale."""
    return dict(os.environ, PATH=tool_path(), LC_ALL="C", **extra)


def tool(name):
    return shutil.which(name, path=tool_path())


def pytest_configure(config):
    config.addinivalue_line("markers", "needs_tools: requires bcftools, samtools and gawk on PATH")


# ---------------------------------------------------------------- fixtures: inputs
@pytest.fixture(scope="session")
def fixtures(tmp_path_factory):
    for t in ("bcftools", "samtools", "bgzip", "tabix"):
        if not tool(t):
            pytest.skip(f"{t} not on PATH")
    out = str(tmp_path_factory.mktemp("fixtures"))
    seqs = mf.build(out)
    return dict(dir=out, seqs=seqs, ref=os.path.join(out, "ref", "test.fa"),
                vcf=os.path.join(out, "vcf"), bam=os.path.join(out, "bam"),
                tumor_bam=os.path.join(out, "bam", "tumor.bam"), normal_bam=os.path.join(out, "bam", "normal.bam"),
                targets=os.path.join(out, "targets.bed"), bseq_header=os.path.join(out, "bseq_header"))


def _dkfz_source():
    for cand in (os.environ.get("SB_DKFZ_DIR"),
                 os.path.join(os.environ.get("SB_DB", "/nonexistent"), "DKFZBiasFilter"),
                 os.path.join(HERE, "_deps", "DKFZBiasFilter")):
        if cand and os.path.isfile(os.path.join(cand, "scripts", "biasFilter.py")):
            return cand
    return None


@pytest.fixture(scope="session")
def mini_db(fixtures, tmp_path_factory):
    """A database directory with only what consensus_scenarios.sh reads: bseq_header and DKFZBiasFilter."""
    db = str(tmp_path_factory.mktemp("db"))
    shutil.copy(fixtures["bseq_header"], os.path.join(db, "bseq_header"))
    src = _dkfz_source()
    if src is None:
        return dict(dir=db, s6=False, reason="DKFZBiasFilter not found (set SB_DKFZ_DIR or SB_DB)")
    dst = os.path.join(db, "DKFZBiasFilter")
    shutil.copytree(src, dst, symlinks=True, ignore=shutil.ignore_patterns(".git"))
    # The public commit a486939 has one Python 2 print statement; patch it so the module imports under Python 3.
    rbf = os.path.join(dst, "scripts", "readbiasfunctions.py")
    try:
        py_compile.compile(rbf, doraise=True)
    except py_compile.PyCompileError:
        text = open(rbf).read().replace('print "No reference sequence for SNV at: "+chrom+" "+str(pos)',
                                        'print("No reference sequence for SNV at: "+chrom+" "+str(pos))')
        open(rbf, "w").write(text)
        py_compile.compile(rbf, doraise=True)
    # Probe: run the filter once on one site exactly as consensus_scenarios.sh does. The public commit is
    # Python 2 code (pysam qqual, map() indexing); if it cannot run under the interpreter given by
    # SB_BIASFILTER_PYTHON (default: the test interpreter), scenario 6 tests are skipped with the reason.
    python = os.environ.get("SB_BIASFILTER_PYTHON", sys.executable)
    probe = os.path.join(db, "probe")
    os.makedirs(os.path.join(probe, "tmp"))
    with open(os.path.join(probe, "in.vcf"), "w") as fh:
        fh.write(open(fixtures["bseq_header"]).read())
        fh.write("chr1\t100\t.\tT\tA\t.\tPASS\t.\t.\n")
    p = subprocess.run([python, os.path.join(dst, "scripts", "biasFilter.py"), os.path.join(probe, "in.vcf"),
                        fixtures["tumor_bam"], fixtures["ref"], os.path.join(probe, "out.vcf"),
                        f"--tempFolder={probe}/tmp", "--mapq=1", "--baseq=1"], capture_output=True, text=True, env=tool_env())
    if p.returncode != 0:
        tail = (p.stderr.strip().splitlines() or ["no output"])[-1]
        return dict(dir=db, s6=False, python=python,
                    reason=f"bias filter does not run under {python}: {tail} (set SB_BIASFILTER_PYTHON)")
    return dict(dir=db, s6=True, python=python, reason="")


# ---------------------------------------------------------------- running the consensus
class Run:
    """One consensus_scenarios.sh invocation and parsed outputs."""

    def __init__(self, workdir, sample, scenarios):
        self.dir = workdir
        self.sample = sample
        self.scenarios = scenarios

    def path(self, sc, kind):
        return os.path.join(self.dir, f"scenario{sc}", f"{self.sample}_{kind}.vcf")

    def records(self, sc, kind):
        """[(chrom, pos, ref, alt, filter, info_dict)] of a scenario output VCF."""
        out = []
        with open(self.path(sc, kind)) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                info = {}
                for kv in f[7].split(";"):
                    if not kv:
                        continue
                    k, _, v = kv.partition("=")
                    info[k] = v if _ else True
                out.append((f[0], int(f[1]), f[3], f[4], f[6], info))
        return out

    def keys(self, sc, kind):
        return {(r[0], r[1], r[2], r[3]) for r in self.records(sc, kind)}

    def record(self, sc, kind, chrom, pos, ref=None, alt=None):
        hits = [r for r in self.records(sc, kind) if r[0] == chrom and r[1] == pos
                and (ref is None or r[2] == ref) and (alt is None or r[3] == alt)]
        assert len(hits) == 1, f"expected one record at {chrom}:{pos} in scenario {sc} {kind}, found {len(hits)}"
        return hits[0]

    def annotated(self, sc, kind):
        """Rows of the annotated TSV as dicts keyed by header name."""
        path = os.path.join(self.dir, f"scenario{sc}", f"{self.sample}_{kind}_annotated.tsv")
        with open(path) as fh:
            header = fh.readline().lstrip("#").rstrip("\n").split("\t")
            rows = [dict(zip(header, l.rstrip("\n").split("\t"))) for l in fh if l.strip()]
        return header, rows

    def decision(self, sc, kind, chrom, pos, ref=None, alt=None):
        header, rows = self.annotated(sc, kind)
        col = next(h for h in header if h.endswith("_DECISION"))
        hits = [r for r in rows if r["CHROM"] == chrom and int(r["POS"]) == pos
                and (ref is None or r["REF"] == ref) and (alt is None or r["ALT"] == alt)]
        assert len(hits) == 1, f"expected one table row at {chrom}:{pos}, found {len(hits)}"
        return hits[0][col]

    def special(self, name):
        return os.path.join(self.dir, "special_mutations", f"{self.sample}_{name}")


def run_consensus(fixtures, mini_db, workdir, scenarios="all", extra=(), muse=None, sample=mf.SAMPLE, expect_fail=False):
    os.makedirs(workdir, exist_ok=True)
    v = fixtures["vcf"]
    cmd = ["bash", CONSENSUS, "-s", scenarios, *extra, sample,
           os.path.join(v, f"{mf.SAMPLE}.mutect2.vcf"),
           muse or os.path.join(v, f"{mf.SAMPLE}.muse.vcf"),
           os.path.join(v, f"{mf.SAMPLE}.strelka_snv.vcf.gz"),
           os.path.join(v, f"{mf.SAMPLE}.strelka_indel.vcf.gz"),
           os.path.join(v, f"{mf.SAMPLE}.sage.vcf.gz"),
           fixtures["tumor_bam"], mini_db["dir"], fixtures["ref"]]
    env = tool_env(PYTHON=mini_db.get("python", sys.executable))
    p = subprocess.run(cmd, cwd=workdir, env=env, capture_output=True, text=True)
    if expect_fail:
        return p
    assert p.returncode == 0, f"consensus_scenarios.sh failed\n--- stdout\n{p.stdout}\n--- stderr\n{p.stderr}"
    run = Run(workdir, sample, scenarios)
    run.stdout, run.stderr = p.stdout, p.stderr
    return run


def _scenarios(mini_db):
    return "all" if mini_db["s6"] else "1,2,3,4,5"


@pytest.fixture(scope="session")
def run_default(fixtures, mini_db, tmp_path_factory):
    return run_consensus(fixtures, mini_db, str(tmp_path_factory.mktemp("run_default")), _scenarios(mini_db))


@pytest.fixture(scope="session")
def run_pon_veto(fixtures, mini_db, tmp_path_factory):
    if not mini_db["s6"]:
        pytest.skip(mini_db["reason"])
    return run_consensus(fixtures, mini_db, str(tmp_path_factory.mktemp("run_pon")), "6", extra=["-p"])


@pytest.fixture(scope="session")
def run_targets(fixtures, mini_db, tmp_path_factory):
    return run_consensus(fixtures, mini_db, str(tmp_path_factory.mktemp("run_targets")), "1,5",
                         extra=["-b", fixtures["targets"]])


@pytest.fixture(scope="session")
def run_contigs(fixtures, mini_db, tmp_path_factory):
    return run_consensus(fixtures, mini_db, str(tmp_path_factory.mktemp("run_contigs")), "1,5",
                         extra=["-c", "chr1,chr2,HPV16"])


@pytest.fixture(scope="session")
def run_empty_muse(fixtures, mini_db, tmp_path_factory):
    return run_consensus(fixtures, mini_db, str(tmp_path_factory.mktemp("run_empty")), "1,5",
                         muse=os.path.join(fixtures["vcf"], "empty.muse.vcf"))


@pytest.fixture(scope="session")
def require_s6(mini_db):
    if not mini_db["s6"]:
        pytest.skip(mini_db["reason"])
