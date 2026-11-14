"""
Proof that the other tests ran against the production tool versions, not merely against whatever is on
PATH. Active only with SB_PRODUCTION_ENVS=1 (the CI job "production-envs" in .github/workflows/ci.yml),
which points SB_TOOLS_BIN, SB_BIASFILTER_PYTHON and SB_VAF_PYTHON at environments created from the
pipeline's own yml/POST.yml (consensus step) and yml/alleleCounter.yml (VAF step). The versions asserted
here are read from those files, so changing a pin in yml/ changes what this test demands.
"""
import glob
import os
import subprocess

import pytest

from conftest import REPO, TOOLS_BIN, VAF_PYTHON, tool, tool_env

pytestmark = pytest.mark.skipif(not os.environ.get("SB_PRODUCTION_ENVS"),
                                reason="SB_PRODUCTION_ENVS not set: tool versions not asserted")

# Python 2 prints `python --version` on stderr, so the version is asked for explicitly.
PY_VERSION = "import sys; print('%d.%d.%d' % sys.version_info[:3])"


def pins(yml):
    yaml = pytest.importorskip("yaml")
    with open(os.path.join(REPO, "yml", yml)) as fh:
        deps = yaml.safe_load(fh)["dependencies"]
    return dict(d.split("=", 1) for d in deps if isinstance(d, str) and "=" in d)


def matches(version, pin):
    """'3.12.7' matches the pins '3', '3.12' and '3.12.7'."""
    return version == pin or version.startswith(pin + ".")


def run(cmd):
    """stdout of a command. stderr is deliberately kept out of the answer: interpreters and tools write
    warnings there (a free-threaded Python re-enabling the GIL to import pysam, for instance)."""
    p = subprocess.run(cmd, capture_output=True, text=True, env=tool_env())
    assert p.returncode == 0, p.stderr
    return p.stdout.strip()


def installed(package):
    """Versions of a conda package installed in the environments listed in SB_TOOLS_BIN (from conda-meta)."""
    found = set()
    for bin_dir in filter(None, TOOLS_BIN.split(":")):
        for meta in glob.glob(os.path.join(bin_dir, "..", "conda-meta", f"{package}-*.json")):
            found.add(os.path.basename(meta)[len(package) + 1:].rsplit("-", 1)[0])
    return found


def test_bcftools_is_the_post_environment_version():
    want = pins("POST.yml")["bcftools"]
    out = run([tool("bcftools"), "--version"])
    assert out.splitlines()[0] == f"bcftools {want}", out


def test_bias_filter_runs_under_the_post_environment_python(mini_db):
    assert mini_db["s6"], mini_db["reason"]
    want = pins("POST.yml")
    version = run([mini_db["python"], "-c", PY_VERSION])
    assert matches(version, want["python"]), version
    assert run([mini_db["python"], "-c", "import pysam; print(pysam.__version__)"]) == want["pysam"]


def test_vaf_runs_under_the_allelecounter_environment_python():
    want = pins("alleleCounter.yml")
    # Without a python pin the solver is free to pick any interpreter (micromamba 2.10 picks a
    # free-threaded Python 3.14), so the environment file must say which Python the VAF step uses.
    assert "python" in want, "yml/alleleCounter.yml does not pin python"
    version = run([VAF_PYTHON, "-c", PY_VERSION])
    assert matches(version, want["python"]), version
    assert run([VAF_PYTHON, "-c", "import pysam; print(pysam.__version__)"]) == want["pysam"]


def test_snv_counts_come_from_the_real_allelecounter():
    # Without alleleCounter on SB_TOOLS_BIN, test_vaf.py falls back to a pysam re-implementation; in the
    # production job the real binary, at the pinned version, must be the one counting.
    want = pins("alleleCounter.yml")["cancerit-allelecount"]
    binary = tool("alleleCounter")
    assert binary, "alleleCounter not found on SB_TOOLS_BIN/PATH"
    assert installed("cancerit-allelecount") == {want}, installed("cancerit-allelecount")
    assert os.path.realpath(os.path.dirname(binary)) in {os.path.realpath(d) for d in TOOLS_BIN.split(":")}, binary
