"""h11: skip the code jobs on prose-only pull requests.

`docs/increments/h11-ci-path-filter.md` §6 specifies these tests:

- T1, T2: the classifier's pure functions in `tools/ci_changes.py`.
- T3: `main` end to end on a fixture git repository.
- T4: the shape of `.github/workflows/main.yaml`, read as text (PyYAML is not
  installed, and not added for this). T4 also runs the `CI result` step's shell
  against the rule of §2.3, since that step is the one required check.
- T5: the suite's own reads. The hook lives in `tests/python/conftest.py`; the
  pure part, `prose_reads`, is unit-tested here and the wiring is proved by a
  planted pytest run in a subprocess.

The tool is loaded lazily through `harness_fixtures.Tool`, so while it is
missing every test that touches it fails naming `tools/ci_changes.py`.
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import sys
from collections.abc import Callable
from pathlib import Path

import pytest

from harness_fixtures import REAL, Tool, clean_env, git

ci = Tool("ci_changes")

WORKFLOW = REAL / ".github" / "workflows" / "main.yaml"
CONFTEST = REAL / "tests" / "python" / "conftest.py"

# ---------------------------------------------------------------------------
# T1: is_prose
# ---------------------------------------------------------------------------

PROSE = (
    "ROADMAP.md",
    "INSTALL.md",
    "testing.md",
    "docs/increments/23c-x.md",
    "docs/retrospectives/2026-10-04-x.md",
    "docs/GLOSSARY.md",
)

CODE = (
    "docs/increments/h4-probes/h4_probe.py",
    "docs/benchmarks/2026-09-26/quarter.geojson",
    ".claude/agents/tester.md",
    ".claude/briefs/common.md",
    ".github/pull_request_template.md",
    ".github/workflows/main.yaml",
    "lib/detria/README.md",
    "tests/python/notes.md",
    "src_python/tin_engine/cli.py",
    "LICENSE",
    "pyproject.toml",
    "CMakeLists.txt",
    ".gitignore",
    "tools/ci_changes.py",
    # Case-sensitive: neither the suffix nor the directory folds case.
    "docs/x.MD",
    "Docs/x.md",
)

#: The Markdown files the build or the suites read (§2.1's table). Every one
#: must be in NOT_PROSE, so the table and the constant cannot drift apart.
READ_BY_CODE = (
    "README.md",
    "CLAUDE.md",
    "docs/increments/README.md",
    "docs/PRINCIPLES.md",
)


@pytest.mark.parametrize("path", PROSE)
def test_t1_prose_paths_are_prose(path: str) -> None:
    assert ci.is_prose(path) is True


@pytest.mark.parametrize("path", CODE)
def test_t1_code_paths_are_not_prose(path: str) -> None:
    assert ci.is_prose(path) is False


@pytest.mark.parametrize("path", READ_BY_CODE)
def test_t1_markdown_read_by_code_is_in_not_prose(path: str) -> None:
    assert path in ci.NOT_PROSE


def test_t1_every_not_prose_entry_is_code() -> None:
    assert ci.NOT_PROSE, "NOT_PROSE is empty"
    for path in ci.NOT_PROSE:
        assert ci.is_prose(path) is False, path


# ---------------------------------------------------------------------------
# T2: needs_full_ci
# ---------------------------------------------------------------------------


def test_t2_empty_list_needs_full_ci() -> None:
    assert ci.needs_full_ci([]) is True


def test_t2_all_prose_does_not_need_full_ci() -> None:
    assert ci.needs_full_ci(list(PROSE)) is False


@pytest.mark.parametrize("code_path", ["src/a.cpp", "README.md", ".gitignore"])
def test_t2_one_code_path_among_prose_needs_full_ci(code_path: str) -> None:
    assert ci.needs_full_ci([*PROSE[:3], code_path, *PROSE[3:]]) is True


# ---------------------------------------------------------------------------
# T3: main end to end on a fixture repository
# ---------------------------------------------------------------------------


def _commit(repo: Path, message: str) -> str:
    git(repo, "add", "-A")
    git(repo, "commit", "-q", "-m", message)
    return git(repo, "rev-parse", "HEAD").strip()


@pytest.fixture
def repo(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """A repository with one code file and one prose file, committed.

    `main` takes only the two SHAs, so it diffs the repository it runs in: the
    test changes into the fixture, and strips GIT_* variables that would
    redirect git elsewhere.
    """
    root = tmp_path / "repo"
    root.mkdir()
    git(root, "init", "-q", "-b", "master")
    git(root, "config", "user.email", "t@example.invalid")
    git(root, "config", "user.name", "T")
    (root / "src").mkdir()
    (root / "src" / "a.cpp").write_text("int main() { return 0; }\n" * 20)
    (root / "docs").mkdir()
    (root / "docs" / "notes.md").write_text("notes\n")
    (root / "ROADMAP.md").write_text("# roadmap\n")
    _commit(root, "root")
    for name in [k for k in os.environ if k.startswith("GIT_")]:
        monkeypatch.delenv(name)
    monkeypatch.chdir(root)
    return root


def _head(repo: Path) -> str:
    return git(repo, "rev-parse", "HEAD").strip()


def _run_main(capsys: pytest.CaptureFixture[str], base: str, head: str) -> tuple[int, str, str]:
    code = ci.main([base, head])
    out, err = capsys.readouterr()
    return code, out, err


def test_t3_prose_only_commit_gives_code_false(
    repo: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    base = _head(repo)
    (repo / "ROADMAP.md").write_text("# roadmap\n\nnew row\n")
    (repo / "docs" / "increments").mkdir()
    (repo / "docs" / "increments" / "h99.md").write_text("design\n")
    head = _commit(repo, "prose")
    assert _run_main(capsys, base, head)[:2] == (0, "code=false\n")


def test_t3_code_commit_gives_code_true(repo: Path, capsys: pytest.CaptureFixture[str]) -> None:
    base = _head(repo)
    (repo / "src" / "a.cpp").write_text("int main() { return 1; }\n")
    (repo / "ROADMAP.md").write_text("# roadmap\n\nrow\n")
    head = _commit(repo, "code and prose")
    assert _run_main(capsys, base, head)[:2] == (0, "code=true\n")


def test_t3_markdown_read_by_code_gives_code_true(
    repo: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    base = _head(repo)
    (repo / "docs" / "increments").mkdir()
    (repo / "docs" / "increments" / "README.md").write_text("protocol\n")
    head = _commit(repo, "protocol text the session-state suite reads")
    assert _run_main(capsys, base, head)[:2] == (0, "code=true\n")


def test_t3_rename_from_code_to_prose_path_gives_code_true(
    repo: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """Red without `--no-renames`: rename detection lists only `docs/a.md`."""
    base = _head(repo)
    git(repo, "mv", "src/a.cpp", "docs/a.md")
    head = _commit(repo, "move a source file to a prose path")
    # Precondition: git itself sees this as a rename, so the probe can fail.
    status = git(repo, "diff", "--name-status", "-M", base, head)
    assert status.startswith("R"), status
    assert _run_main(capsys, base, head)[:2] == (0, "code=true\n")


def test_t3_changed_paths_lists_both_sides_of_a_rename(repo: Path) -> None:
    base = _head(repo)
    git(repo, "mv", "src/a.cpp", "docs/a.md")
    head = _commit(repo, "rename")
    assert set(ci.changed_paths(base, head, repo)) == {"src/a.cpp", "docs/a.md"}


def test_t3_deleting_a_source_file_alone_gives_code_true(
    repo: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    base = _head(repo)
    git(repo, "rm", "-q", "src/a.cpp")
    head = _commit(repo, "delete")
    assert _run_main(capsys, base, head)[:2] == (0, "code=true\n")


def test_t3_unusual_prose_file_name_is_read_unquoted(
    repo: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """Red without `-z`: git quotes the name, and `"..."` does not end in `.md`."""
    base = _head(repo)
    (repo / "docs" / "naïve notes.md").write_text("prose\n")
    head = _commit(repo, "unusual name")
    assert _run_main(capsys, base, head)[:2] == (0, "code=false\n")


@pytest.mark.parametrize("base_kind", ["unknown SHA", "empty base", "base equals head"])
def test_t3_unusable_base_falls_back_to_a_full_run(
    repo: Path, capsys: pytest.CaptureFixture[str], base_kind: str
) -> None:
    (repo / "ROADMAP.md").write_text("# roadmap\n\nrow\n")
    head = _commit(repo, "prose")
    base = {
        "unknown SHA": "0123456789abcdef0123456789abcdef01234567",
        "empty base": "",
        "base equals head": head,
    }[base_kind]
    code, out, err = _run_main(capsys, base, head)
    assert (code, out) == (0, "code=true\n")
    assert err.strip(), f"{base_kind}: no sentence on stderr saying why"


# ---------------------------------------------------------------------------
# T4: the workflow's shape
# ---------------------------------------------------------------------------

#: Jobs that do not gate on `changes`: the classifier itself, governance (runs
#: on every event), and the aggregate.
UNGATED = ("changes", "governance", "result")
#: Jobs left out of `CI result`'s needs on purpose (§9, ruling 1).
EXEMPT_FROM_RESULT = ("tsan",)
#: Jobs on master before h11; the parser must find every one, or it is parsing
#: a layout it does not understand.
KNOWN_JOBS = {"cpp", "sanitizers", "tsan", "governance", "python"}

GATE = "if: needs.changes.outputs.code != 'false'"


def _top_level_block(text: str, key: str) -> list[str]:
    """The lines under a top-level `key:` up to the next top-level key."""
    lines = text.splitlines()
    starts = [i for i, line in enumerate(lines) if line.rstrip() == f"{key}:"]
    assert len(starts) == 1, f"main.yaml: expected one top-level '{key}:' line"
    block = []
    for line in lines[starts[0] + 1 :]:
        if line and not line[0].isspace() and not line.startswith("#"):
            break
        block.append(line)
    return block


def _jobs(text: str) -> dict[str, list[str]]:
    """Each job's lines, keyed by job id (two-space indent under `jobs:`)."""
    jobs: dict[str, list[str]] = {}
    current: str | None = None
    for line in _top_level_block(text, "jobs"):
        match = re.fullmatch(r"  ([A-Za-z_][\w-]*):\s*", line)
        if match:
            current = match.group(1)
            jobs[current] = []
        elif current is not None:
            jobs[current].append(line)
    missing = KNOWN_JOBS - jobs.keys()
    assert not missing, f"main.yaml: jobs not found by the parser: {sorted(missing)}"
    return jobs


def _job_keys(lines: list[str]) -> dict[str, str]:
    """A job's own keys (four-space indent) and their inline values."""
    keys = {}
    for line in lines:
        match = re.fullmatch(r"    ([\w-]+):\s*(.*?)\s*", line)
        if match:
            keys[match.group(1)] = match.group(2)
    return keys


def _flow_list(value: str) -> set[str]:
    match = re.fullmatch(r"\[(.*)\]", value)
    assert match, f"expected a one-line [a, b] list, got {value!r}"
    return {item.strip() for item in match.group(1).split(",") if item.strip()}


@pytest.fixture(scope="module")
def workflow() -> str:
    return WORKFLOW.read_text()


@pytest.fixture(scope="module")
def jobs(workflow: str) -> dict[str, list[str]]:
    return _jobs(workflow)


def test_t4_triggers_have_no_push(workflow: str) -> None:
    triggers = {
        m.group(1)
        for line in _top_level_block(workflow, "on")
        if (m := re.fullmatch(r"  ([\w-]+):.*", line))
    }
    assert {"pull_request", "merge_group"} <= triggers, triggers
    assert "workflow_dispatch" in triggers
    assert "push" not in triggers


def test_t4_changes_and_result_jobs_exist_under_their_names(
    jobs: dict[str, list[str]],
) -> None:
    assert _job_keys(jobs.get("changes", [])).get("name") == "Changed files"
    # `CI result` is the one required check (§2.3); its name is what branch
    # protection matches, so it is pinned exactly.
    assert _job_keys(jobs.get("result", [])).get("name") == "CI result"


def test_t4_changes_job_runs_the_classifier_with_full_history(
    jobs: dict[str, list[str]],
) -> None:
    body = "\n".join(jobs.get("changes", []))
    assert "fetch-depth: 0" in body
    assert "tools/ci_changes.py" in body
    assert '>> "$GITHUB_OUTPUT"' in body
    assert re.search(r"^    outputs:\s*$", body, re.M), "changes has no outputs:"
    assert re.search(r"^      code: \$\{\{ steps\.[\w-]+\.outputs\.code \}\}", body, re.M)


def test_t4_every_other_job_is_gated_on_changes(jobs: dict[str, list[str]]) -> None:
    gated = [job for job in jobs if job not in UNGATED]
    assert {"cpp", "sanitizers", "tsan", "python"} <= set(gated)
    for job in gated:
        keys = _job_keys(jobs[job])
        assert keys.get("needs") == "changes", f"{job}: needs is {keys.get('needs')!r}"
        assert f"if: {keys.get('if')}" == GATE, f"{job}: if is {keys.get('if')!r}"


def test_t4_governance_runs_on_every_event(jobs: dict[str, list[str]]) -> None:
    keys = _job_keys(jobs["governance"])
    assert "needs" not in keys and "if" not in keys, keys


def test_t4_result_always_runs_and_needs_every_gating_job(
    jobs: dict[str, list[str]],
) -> None:
    keys = _job_keys(jobs.get("result", []))
    assert keys.get("if") == "always()", keys
    expected = set(jobs) - {"result"} - set(EXEMPT_FROM_RESULT)
    assert _flow_list(keys.get("needs", "")) == expected, (
        "a job must be in CI result's needs or in EXEMPT_FROM_RESULT"
    )


def test_t4_unchecked_cpp_leg_does_not_block(jobs: dict[str, list[str]]) -> None:
    keys = _job_keys(jobs["cpp"])
    assert keys.get("continue-on-error") == "${{ matrix.hardening == 'OFF' }}", keys


# T4, continued: the `CI result` step's shell, run against §2.3's rule.


def _result_script(lines: list[str]) -> str:
    """The `run: |` block of the result job's one step, dedented."""
    starts = [i for i, line in enumerate(lines) if re.fullmatch(r"\s+run: \|.*", line)]
    assert len(starts) == 1, "result: expected one 'run: |' block"
    indent = len(lines[starts[0]]) - len(lines[starts[0]].lstrip()) + 2
    script = []
    for line in lines[starts[0] + 1 :]:
        if line.strip() and len(line) - len(line.lstrip()) < indent:
            break
        script.append(line[indent:])
    text = "\n".join(script).strip()
    assert text, "result: the run block is empty"
    assert "${{" not in text, "result: inputs come through env, not ${{ }} in run"
    return text + "\n"


S, K, F, C = "success", "skipped", "failure", "cancelled"


def _verdict(script: str, code: str, changes: str, governance: str, others: str) -> tuple[int, str]:
    bash = shutil.which("bash")
    assert bash is not None
    env = {
        **clean_env(),
        "CODE": code,
        "CHANGES": changes,
        "GOVERNANCE": governance,
        "OTHERS": others,
    }
    # GitHub runs `run:` with `bash --noprofile --norc -eo pipefail {0}`.
    done = subprocess.run(
        [bash, "--noprofile", "--norc", "-eo", "pipefail", "-c", script],
        env=env,
        capture_output=True,
        text=True,
        timeout=30,
    )
    return done.returncode, done.stdout + done.stderr


@pytest.mark.parametrize(
    ("code", "changes", "governance", "others"),
    [
        ("true", S, S, f"{S} {S} {S}"),
        ("false", S, S, f"{K} {K} {K}"),
        ("false", S, S, f"{S} {K} {S}"),
    ],
    ids=["code change, all green", "prose, all skipped", "prose, mixed"],
)
def test_t4_result_passes(
    jobs: dict[str, list[str]], code: str, changes: str, governance: str, others: str
) -> None:
    result_script = _result_script(jobs.get("result", []))
    status, output = _verdict(result_script, code, changes, governance, others)
    assert status == 0, output


@pytest.mark.parametrize(
    ("code", "changes", "governance", "others", "named"),
    [
        ("true", S, S, f"{S} {K} {S}", K),
        ("", S, S, f"{K} {K} {K}", K),
        ("garbage", S, S, f"{K} {K} {K}", K),
        ("false", S, K, f"{K} {K} {K}", K),
        ("false", S, F, f"{K} {K} {K}", F),
        ("", F, S, f"{K} {K} {K}", F),
        ("true", S, S, f"{S} {F} {S}", F),
        ("false", S, S, f"{K} {K} {F}", F),
        ("true", S, S, f"{S} {S} {C}", C),
        ("false", C, S, f"{K} {K} {K}", C),
    ],
    ids=[
        "code change, a job skipped",
        "empty code output, jobs skipped",
        "unexpected code output, jobs skipped",
        "governance skipped",
        "governance failed",
        "changes job crashed",
        "code change, a job failed",
        "prose, a job failed",
        "a job cancelled",
        "changes cancelled",
    ],
)
def test_t4_result_fails_and_says_how(
    jobs: dict[str, list[str]],
    code: str,
    changes: str,
    governance: str,
    others: str,
    named: str,
) -> None:
    result_script = _result_script(jobs.get("result", []))
    status, output = _verdict(result_script, code, changes, governance, others)
    assert status == 1, output
    assert named in output, f"the failure does not say which job ended {named!r}"


# ---------------------------------------------------------------------------
# T5: the suite's own reads
# ---------------------------------------------------------------------------


def test_t5_prose_reads_names_prose_inside_the_root(tmp_path: Path) -> None:
    root = tmp_path / "root"
    paths = [
        str(root / "ROADMAP.md"),
        str(root / "docs" / "increments" / "h99.md"),
        str(root / "README.md"),  # in NOT_PROSE
        str(root / "docs" / "increments" / "README.md"),  # in NOT_PROSE
        str(root / "tools" / "ci_changes.py"),
        str(root / ".claude" / "agents" / "tester.md"),
        str(tmp_path / "elsewhere" / "ROADMAP.md"),  # outside the root
        str(tmp_path / "ROADMAP.md"),  # outside the root
    ]
    assert set(ci.prose_reads(paths, root)) == {"ROADMAP.md", "docs/increments/h99.md"}


def test_t5_prose_reads_of_nothing_is_empty(tmp_path: Path) -> None:
    assert list(ci.prose_reads([], tmp_path)) == []


PLANTED = """\
from pathlib import Path

REAL = Path({real!r})


def test_planted_read() -> None:
    assert (REAL / {target!r}).read_text()
"""


@pytest.fixture
def planted(tmp_path: Path) -> Callable[[str], subprocess.CompletedProcess[str]]:
    """Run one planted test in a subprocess, with the suite's hook loaded.

    The planted file sits outside `tests/python`, so the conftest is loaded as a
    plugin (`-p conftest`, `tests/python` on the path); `addopts` is cleared so
    coverage settings do not apply to the scratch run.
    """

    def run(target: str) -> subprocess.CompletedProcess[str]:
        assert CONFTEST.exists(), "tests/python/conftest.py (the T5 hook) does not exist"
        test = tmp_path / f"test_planted_{target.replace('/', '_').replace('.', '_')}.py"
        test.write_text(PLANTED.format(real=str(REAL), target=target))
        env = {**clean_env(), "PYTHONPATH": str(CONFTEST.parent)}
        return subprocess.run(
            [
                sys.executable,
                "-m",
                "pytest",
                "-q",
                "-p",
                "conftest",
                "-o",
                "addopts=",
                "--rootdir",
                str(tmp_path),
                str(test),
            ],
            cwd=tmp_path,
            env=env,
            capture_output=True,
            text=True,
            timeout=120,
        )

    return run


def test_t5_a_suite_that_reads_prose_fails_the_session(
    planted: Callable[[str], subprocess.CompletedProcess[str]],
) -> None:
    done = planted("ROADMAP.md")
    output = done.stdout + done.stderr
    assert done.returncode != 0, output
    assert "ROADMAP.md" in output, output
    assert "NOT_PROSE" in output, output


def test_t5_a_suite_that_reads_only_code_passes(
    planted: Callable[[str], subprocess.CompletedProcess[str]],
) -> None:
    """The control: the planted run itself works, so the failure above is the hook's."""
    done = planted("pyproject.toml")
    assert done.returncode == 0, done.stdout + done.stderr
