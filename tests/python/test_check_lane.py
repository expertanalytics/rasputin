"""`tools/check_lane.py` and the two CI steps that hold the lanes and the size ceiling.

Spec: `docs/increments/h34-lean-process.md` §3 (*Enforced in CI*) and §4,
tests 11 and 13. A branch that changes a full-lane path (the meshing core,
elevation sampling, hydrology, catchments, map projection) must add an
increment file, `docs/increments/<name>.md` one level down, `README.md`
excluded; otherwise the tool names each such path and exits 1. A git failure
exits 2. Paths come from `git diff --name-only --no-renames <base>...<head>`,
so a move out of a full-lane directory still names the path it left.

The tool runs as `python3 tools/check_lane.py <base> <head>` from a copy
committed in the scratch repository's root commit, so the copy is never part
of the diff it reads. Test 13 reads `.github/workflows/main.yaml` as text.
"""

from __future__ import annotations

import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from harness_fixtures import REAL, clean_env, git

pytestmark = pytest.mark.harness

SCRIPT = "tools/check_lane.py"

#: §3's `FULL_LANE`: one changed file under or at each entry.
FULL_LANE_FILES = (
    "include/terrain/cdt/x.hpp",
    "include/terrain/predicates/x.hpp",
    "include/terrain/refinement/x.hpp",
    "include/terrain/noding/x.hpp",
    "include/terrain/core/x.hpp",
    "include/terrain/mesh/x.hpp",
    "src/cdt/x.cpp",
    "include/terrain/raster/x.hpp",
    "src_python/tin_engine/elevation.py",
    "src_python/tin_engine/raster.py",
    "src_python/tin_engine/dem_input.py",
    "src_python/tin_engine/mosaic.py",
    "src_python/tin_engine/target_grid.py",
    "include/terrain/hydrology/x.hpp",
    "src_python/tin_engine/catchment.py",
    "src_python/tin_engine/catchment_core.py",
    "src_python/tin_engine/network.py",
    "src_python/tin_engine/sub_catchments.py",
    "src_python/tin_engine/crs.py",
    "src_python/tin_engine/geodesy.py",
)

#: Fast-lane paths: none starts with a `FULL_LANE` entry. `src_python/` shares
#: `src`'s first letters but not `src/`.
FAST_LANE_FILES = (
    "src_python/tin_engine/cli.py",
    "src_python/tin_engine/catchment_batch.py",
    "include/terrain/land_cover/x.hpp",
    "tools/duties.py",
    "tests/python/test_x.py",
    "docs/x.md",
    "ROADMAP.md",
)

INCREMENT = "docs/increments/99-x.md"


def commit_files(repo: Path, files: dict[str, str | None], message: str) -> None:
    """Write (or, for None, delete) each file and commit everything."""
    for relative, text in files.items():
        path = repo / relative
        if text is None:
            path.unlink()
        else:
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(text)
    git(repo, "add", "-A")
    git(repo, "commit", "-q", "-m", message)


@pytest.fixture
def repo(tmp_path: Path) -> Path:
    """A repository on `master` whose root commit holds the script's copy and
    one file in each of the core's directories; `branch` starts there."""
    root = tmp_path.resolve() / "repo"
    root.mkdir()
    git(root, "init", "-q", "-b", "master")
    git(root, "config", "user.email", "t@example.invalid")
    git(root, "config", "user.name", "T")
    git(root, "config", "commit.gpgsign", "false")
    source = REAL / SCRIPT
    if source.exists():
        (root / SCRIPT).parent.mkdir(parents=True)
        shutil.copy2(source, root / SCRIPT)
    commit_files(root, {"include/terrain/core/a.hpp": "int a;\n", "README.md": "# r\n"}, "base")
    git(root, "checkout", "-q", "-b", "branch")
    return root


def run_lane(
    repo: Path, base: str = "master", head: str = "branch"
) -> subprocess.CompletedProcess[str]:
    script = repo / SCRIPT
    if not script.exists():
        pytest.fail(f"{SCRIPT} is missing from the checkout")
    return subprocess.run(
        [sys.executable, str(script), base, head],
        capture_output=True,
        text=True,
        cwd=repo,
        env=clean_env(),
        timeout=60,
        check=False,
    )


def on_branch(repo: Path, *paths: str) -> subprocess.CompletedProcess[str]:
    commit_files(repo, {path: f"// {path}\n" for path in paths}, "change")
    return run_lane(repo)


def refused_naming(result: subprocess.CompletedProcess[str], *paths: str) -> str:
    """Exit 1, and each of `paths` named in the output."""
    said = result.stdout + result.stderr
    assert result.returncode == 1, f"exit {result.returncode}, not 1; output: {said}"
    for path in paths:
        assert path in said, f"{path} is not named; output: {said}"
    return said


def passed(result: subprocess.CompletedProcess[str]) -> None:
    said = result.stdout + result.stderr
    assert result.returncode == 0, f"exit {result.returncode}, not 0; output: {said}"


# ---------------------------------------------------------------- test 11


def test_h34_11_a_noding_change_alone_is_refused_and_named(repo: Path) -> None:
    refused_naming(on_branch(repo, "include/terrain/noding/x.hpp"), "include/terrain/noding/x.hpp")


def test_h34_11_a_noding_change_with_an_increment_file_passes(repo: Path) -> None:
    passed(on_branch(repo, "include/terrain/noding/x.hpp", INCREMENT))


def test_h34_11_the_increments_readme_is_not_an_increment_file(repo: Path) -> None:
    refused_naming(
        on_branch(repo, "include/terrain/noding/x.hpp", "docs/increments/README.md"),
        "include/terrain/noding/x.hpp",
    )


def test_h34_11_an_increment_file_two_levels_down_does_not_count(repo: Path) -> None:
    """§3: `docs/increments/<name>.md`, one level."""
    refused_naming(
        on_branch(repo, "include/terrain/noding/x.hpp", "docs/increments/old/99-x.md"),
        "include/terrain/noding/x.hpp",
    )


def test_h34_11_only_the_full_lane_paths_are_named(repo: Path) -> None:
    said = refused_naming(
        on_branch(repo, "include/terrain/noding/x.hpp", "src_python/tin_engine/cli.py"),
        "include/terrain/noding/x.hpp",
    )
    assert "src_python/tin_engine/cli.py" not in said, said


@pytest.mark.parametrize("path", FULL_LANE_FILES)
def test_h34_11_each_full_lane_path_alone_is_refused(repo: Path, path: str) -> None:
    refused_naming(on_branch(repo, path), path)


@pytest.mark.parametrize("path", FAST_LANE_FILES)
def test_h34_11_a_fast_lane_path_alone_passes(repo: Path, path: str) -> None:
    passed(on_branch(repo, path))


def test_h34_11_an_unknown_base_is_exit_2(repo: Path) -> None:
    commit_files(repo, {"src_python/tin_engine/cli.py": "x = 1\n"}, "change")
    result = run_lane(repo, base="no-such-base")
    assert result.returncode == 2, (
        f"exit {result.returncode}; output: {result.stdout}{result.stderr}"
    )
    assert "no-such-base" in result.stderr, result.stderr


def test_h34_11_a_move_out_of_a_full_lane_directory_is_refused(repo: Path) -> None:
    """`--no-renames`: the move lists both paths, so the one it left is named."""
    (repo / "src_python").mkdir()
    git(repo, "mv", "include/terrain/core/a.hpp", "src_python/a.hpp")
    git(repo, "commit", "-q", "-m", "move")
    status = git(repo, "diff", "-M", "--name-status", "master", "branch")
    assert status.startswith("R"), f"git did not see a rename:\n{status}"
    refused_naming(run_lane(repo), "include/terrain/core/a.hpp")


def test_h34_11_the_diff_is_from_the_merge_base(repo: Path) -> None:
    """`<base>...<head>`: a full-lane change master made after the branch
    started is not the branch's."""
    commit_files(repo, {"src_python/tin_engine/cli.py": "x = 1\n"}, "fast change")
    git(repo, "checkout", "-q", "master")
    commit_files(repo, {"include/terrain/noding/y.hpp": "int y;\n"}, "master moves on")
    git(repo, "checkout", "-q", "branch")
    passed(run_lane(repo))


# ---------------------------------------------------------------- test 13, the workflow

WORKFLOW = REAL / ".github" / "workflows" / "main.yaml"
BASE = "${{ github.event.pull_request.base.sha }}"
HEAD = "${{ github.sha }}"
ON_PULL_REQUEST = "if: github.event_name == 'pull_request'"


def governance_steps() -> list[str]:
    """The governance job's steps, each as its lines joined, whitespace collapsed."""
    lines = WORKFLOW.read_text().splitlines()
    starts = [i for i, line in enumerate(lines) if line.rstrip() == "  governance:"]
    assert len(starts) == 1, "main.yaml: expected one governance job"
    job: list[str] = []
    for line in lines[starts[0] + 1 :]:
        if re.fullmatch(r"  [A-Za-z_][\w-]*:\s*", line):
            break
        job.append(line)
    steps: list[list[str]] = []
    for line in job:
        if re.match(r"      - ", line):
            steps.append([])
        if steps and not line.lstrip().startswith("#"):
            steps[-1].append(line)
    return [" ".join(" ".join(step).split()) for step in steps]


def step_running(tool: str) -> str:
    found = [step for step in governance_steps() if tool in step]
    assert len(found) == 1, f"{len(found)} governance steps run {tool}"
    return found[0]


@pytest.mark.parametrize("tool", ["tools/check_lane.py", "tools/count_loc.py"])
def test_h34_13_the_governance_job_runs_the_step_on_pull_requests_only(tool: str) -> None:
    step = step_running(tool)
    assert ON_PULL_REQUEST in step, step


@pytest.mark.parametrize("tool", ["tools/check_lane.py", "tools/count_loc.py"])
def test_h34_13_the_step_reads_the_pull_requests_base_then_its_head(tool: str) -> None:
    step = step_running(tool)
    run = step[step.index(tool) :]
    assert BASE in run and HEAD in run, step
    assert run.index(BASE) < run.index(HEAD), step
    assert "merge_group" not in step, step


def test_h34_13_the_size_step_holds_the_ceiling_at_700() -> None:
    step = step_running("tools/count_loc.py")
    assert re.search(r"--max[ =]700\b", step), step
