"""`tools/push_grant.py`: the granted forms and the reviewed rule (h20 §4.2, §4.4).

The spec is `docs/increments/h20-away-push-grants.md` §3.2 (the three forms and
the expansion test) and §3.3 (`reviewed(H)`). §4.4 is the increment's
invariant-critical suite: each case builds on `harness_fixtures.Granted`, and
the expected result is a pass (`None`) or the reason the design names.

Pinned here beyond the design's words: `reviewed(repo, head, records)` takes
`records` as the parsed lines of `reviews.jsonl`, oldest first; a short sha
in a reason is at least 7 hex digits; `grants(flag)` returns a tuple.
"""

from __future__ import annotations

import json
import os
import re
import subprocess
from collections.abc import Callable
from pathlib import Path
from typing import Any

import pytest

from harness_fixtures import (
    PRIVATE_REPO,
    Granted,
    PrivateGranted,
    Tool,
    commit_file,
    git,
    grants_used_path,
    jsonl,
    make_repo,
    review_record,
    write_reviews,
)

# h17 §4b: a harness test; CI runs it in the `harness` job, not the product legs.
pytestmark = pytest.mark.harness

push_grant = Tool("push_grant")

H = "0123456789abcdef0123456789abcdef01234567"
SHORT = H[:39]
PUSH = f"git push origin {H}:refs/heads/feat"
CREATE = "gh pr create --head feat --title 't' --body-file f"
MERGE = f"gh pr merge feat --match-head-commit {H} --merge"

NO_REVIEW = "no review recorded for this branch"
CHANGES = "last review round: CHANGES REQUESTED"


# ---------------------------------------------------------------- test 8

FORMS = (
    PUSH,
    f"git -C /abs/dir push origin {H}:refs/heads/feat",
    "gh pr create --head feat",
    CREATE,
    "gh pr create -H feat -B master -t 'a title' -F f -d",
    "gh pr create --head=feat --base=master --title 'a title' --body-file=f --draft",
    "gh pr create --body-file f --base master --title 'a title' --head feat",
    f"gh pr merge feat --match-head-commit {H}",
    MERGE,
    f"gh pr merge --merge --match-head-commit {H} feat",
)


@pytest.mark.parametrize("command", FORMS)
def test_each_granted_form_parses(command: str) -> None:
    assert push_grant.granted_act(command) is not None


# ---------------------------------------------------------------- test 9

NOT_FORMS = (
    f"git push --force origin {H}:refs/heads/feat",
    f"git push -f origin {H}:refs/heads/feat",
    f"git push --force-with-lease origin {H}:refs/heads/feat",
    f"git push --no-verify origin {H}:refs/heads/feat",
    f"git push -u origin {H}:refs/heads/feat",
    "git push origin feat",
    f"git push origin +{H}:refs/heads/feat",
    f"git push origin {H}:feat",
    f"git push origin {H}:refs/heads/main",
    "git push origin HEAD",
    f"git push upstream {H}:refs/heads/feat",
    f"git push origin {H}:refs/heads/feat {H}:refs/heads/other",
    f"git push origin {SHORT}:refs/heads/feat",
    f"gh pr merge feat --match-head-commit {H} --admin",
    f"gh pr merge feat --match-head-commit {H} -d",
    f"gh pr merge feat --match-head-commit {H} --auto",
    f"gh pr merge feat --match-head-commit {H} --squash",
    f"gh pr merge feat --match-head-commit {H} --rebase",
    f"gh pr merge feat --match-head-commit {H} -R x/y",
    f"gh pr merge feat --match-head-commit {SHORT}",
    "gh pr merge feat --merge",
    f"gh pr merge 241 --match-head-commit {H}",
    "gh pr create --head feat --body text",
    "gh pr create --head feat --base other",
    # The expansions: words the shell would change, which shlex.split does not.
    f"git -C ~/other push origin {H}:refs/heads/feat",
    f"git -C =other push origin {H}:refs/heads/feat",
    f"git -C rel/dir push origin {H}:refs/heads/feat",
    f"git -C '$HOME/x' push origin {H}:refs/heads/feat",
    "gh pr create --head 'fea*'",
    "gh pr create --head feat --body-file ~/f",
    "gh pr create --head feat --body-file==ls",
    "gh pr create --head feat --body-file=~/f",
    "gh pr create --head feat --title t",
    'gh pr create --head feat --title "t"',
)


@pytest.mark.parametrize("command", NOT_FORMS)
def test_anything_outside_the_forms_is_not_a_granted_act(command: str) -> None:
    assert push_grant.granted_act(command) is None


# ---------------------------------------------------------------- grants (§3.1)


def test_the_grants_of_a_flag_are_its_list() -> None:
    assert tuple(push_grant.grants({"grants": ["feat", "other"]})) == ("feat", "other")


@pytest.mark.parametrize(
    "flag",
    [{}, {"grants": None}, {"grants": "feat"}, {"grants": ["feat", 1]}, {"grants": {"feat": 1}}],
    ids=["written-before-h20", "null", "a-string", "a-non-string-item", "a-mapping"],
)
def test_a_flag_without_a_list_of_strings_grants_nothing(flag: dict[str, object]) -> None:
    assert tuple(push_grant.grants(flag)) == ()


# ---------------------------------------------------------------- §4.4


@pytest.fixture
def built(tmp_path: Path) -> Granted:
    return Granted(tmp_path / "repo", flag=False)


def judge(built: Granted) -> str | None:
    found: str | None = push_grant.reviewed(built.repo, built.head, built.records())
    return found


def refuses_naming(reason: str | None, commit: str, path: str, review: str) -> None:
    assert reason is not None, "passed; expected a refusal"
    pattern = (
        rf"commit {commit[:7]}[0-9a-f]* changes {re.escape(path)} after the review of {review[:7]}"
    )
    assert re.search(pattern, reason), reason


def test_19_the_fixture_as_built_passes(built: Granted) -> None:
    assert judge(built) is None


def test_20_no_record_is_no_review(built: Granted) -> None:
    write_reviews(built.repo, [])
    assert NO_REVIEW in (judge(built) or "")


def test_20_a_record_of_another_branch_is_no_review(built: Granted) -> None:
    git(built.repo, "checkout", "-q", "-b", "other", "master")
    elsewhere = commit_file(built.repo, "tools/s.py", "s = 1\n", "S on other")
    git(built.repo, "checkout", "-q", "feat")
    write_reviews(built.repo, [review_record(elsewhere)])
    assert NO_REVIEW in (judge(built) or "")


@pytest.mark.parametrize(
    ("first", "second", "expected"),
    [
        (("x", "APPROVED"), ("x", "CHANGES REQUESTED"), CHANGES),
        (("x", "CHANGES REQUESTED"), ("x", "APPROVED"), None),
        (("x", "APPROVED"), ("rec", "CHANGES REQUESTED"), CHANGES),
        (("x", "CHANGES REQUESTED"), ("rec", "APPROVED"), None),
        # The newest is the last line of the file, not the record nearest H.
        (("rec", "APPROVED"), ("x", "CHANGES REQUESTED"), CHANGES),
    ],
    ids=["same-sha-refused", "same-sha-passes", "approved-then-changes", "changes-then-approved",
         "last-line-wins-over-nearest-head"],
)  # fmt: skip
def test_21_the_newest_record_decides(
    built: Granted, first: tuple[str, str], second: tuple[str, str], expected: str | None
) -> None:
    sha = {"x": built.x, "rec": built.rec}
    write_reviews(built.repo, [review_record(sha[first[0]], first[1]),
                               review_record(sha[second[0]], second[1])])  # fmt: skip
    found = judge(built)
    if expected is None:
        assert found is None
    else:
        assert found is not None and expected in found, found


def test_22_approved_with_wordings_and_a_prose_recording_commit_pass(built: Granted) -> None:
    write_reviews(built.repo, [review_record(built.x, "APPROVED WITH WORDINGS")])
    assert judge(built) is None


def test_23_a_code_commit_after_the_recording_commit_is_refused(built: Granted) -> None:
    late = commit_file(built.repo, "tools/x.py", "x = 3\n", "a late code change")
    refuses_naming(judge(built), late, "tools/x.py", built.x)


def test_24_a_recording_commit_that_changes_code_is_refused(tmp_path: Path) -> None:
    built = Granted(tmp_path / "repo", flag=False, recording_code=True)
    refuses_naming(judge(built), built.rec, "tools/x.py", built.x)


@pytest.mark.parametrize("path", ["CLAUDE.md", "README.md", "docs/increments/README.md"])
def test_25_a_post_review_change_to_a_markdown_file_the_suites_read_is_refused(
    built: Granted, path: str
) -> None:
    late = commit_file(built.repo, path, "# changed after the review\n", f"edit {path}")
    refuses_naming(judge(built), late, path, built.x)


def test_26_a_clean_master_merge_passes(built: Granted) -> None:
    built.on_master("src/y.py", "y = 1\n")
    built.merge("refs/remotes/origin/master")
    assert judge(built) is None


def test_27_a_master_merge_with_a_roadmap_conflict_resolved_by_hand_passes(
    built: Granted,
) -> None:
    commit_file(built.repo, "ROADMAP.md", "| feat | done |\n", "ROADMAP row for feat")
    built.on_master("ROADMAP.md", "| other | done |\n")
    built.merge("refs/remotes/origin/master",
                resolve={"ROADMAP.md": "| other | done |\n| feat | done |\n"})  # fmt: skip
    assert judge(built) is None


def test_28_a_master_merge_with_a_code_conflict_resolved_by_hand_is_refused(
    built: Granted,
) -> None:
    built.on_master("tools/x.py", "x = 'master'\n")
    merge = built.merge("refs/remotes/origin/master", resolve={"tools/x.py": "x = 'both'\n"})
    refuses_naming(judge(built), merge, "tools/x.py", built.x)


def test_29_a_clean_merge_of_a_branch_not_on_origin_master_is_refused(built: Granted) -> None:
    git(built.repo, "checkout", "-q", "-b", "side", "master")
    commit_file(built.repo, "src/z.py", "z = 1\n", "side changes code")
    git(built.repo, "checkout", "-q", "feat")
    merge = built.merge("side")
    found = judge(built)
    assert found is not None, "a merge of a non-master branch passed as a master merge"
    assert f"commit {merge[:7]}" in found, found


def test_30_a_review_on_a_merged_branch_off_the_first_parent_chain_is_not_seen(
    tmp_path: Path,
) -> None:
    """feat merges `other` (code commit S) before X; a newer APPROVED at S must not count.

    A full `origin/master..H` walk picks S's record and refuses; the
    first-parent walk picks X's and passes.
    """

    repo = make_repo(tmp_path / "repo")
    git(repo, "update-ref", "refs/remotes/origin/master", "master")
    git(repo, "checkout", "-q", "-b", "other")
    s = commit_file(repo, "tools/s.py", "s = 1\n", "S on other")
    git(repo, "checkout", "-q", "-b", "feat", "master")
    git(repo, "merge", "-q", "--no-ff", "--no-edit", "other")
    x = commit_file(repo, "tools/x.py", "x = 1\n", "X: a code commit")
    commit_file(repo, "docs/increments/h99.md", "review recorded\n", "h99 review recorded")
    head = git(repo, "rev-parse", "feat").strip()
    records = [review_record(x), review_record(s)]
    assert push_grant.reviewed(repo, head, records) is None


def test_31_a_review_of_masters_own_commit_is_not_this_branchs(built: Granted) -> None:
    m = built.on_master("src/m.py", "m = 1\n")
    write_reviews(built.repo, [review_record(m)])
    built.merge("refs/remotes/origin/master")
    assert NO_REVIEW in (judge(built) or "")


def test_32_a_review_line_in_a_commit_message_is_not_read(built: Granted) -> None:
    """Y's message cannot name Y, so the prose commit after Y names it instead."""
    y = commit_file(built.repo, "tools/x.py", "x = 4\n", "Y: a code change")
    commit_file(built.repo, "docs/increments/h99.md", "again\n", f"note\n\nReview: APPROVED {y}")
    refuses_naming(judge(built), y, "tools/x.py", built.x)


# ================================================================ h22 fix 1
#
# docs/increments/h22-harness-fixes.md §4 and §7 tests 1-8: the private forms,
# a branch's own repository (`home`), the walk against the branch's own
# master, and the checks of a private merge. Every case runs on Fixture A
# (`harness_fixtures.PrivateGranted`) in `tmp_path`; `gh` is never run: the
# one `subprocess.run` that `push_grant` makes for it is answered by `Runner`,
# and a `gh` on PATH that fails loudly backs that up.

P = PRIVATE_REPO
PRIVATE_FORMS = (
    (f"git push private {H}:refs/heads/prv", "push", H, None),
    (f"git -C /abs/dir push private {H}:refs/heads/prv", "push", H, "/abs/dir"),
    (f"gh pr create -R {P} --head prv", "pull request", None, None),
    (f"gh pr create --repo {P} --head prv --base master --title 't' --body-file f --draft",
     "pull request", None, None),
    (f"gh pr create --repo={P} --head prv", "pull request", None, None),
    (f"gh pr merge prv -R {P} --merge --match-head-commit {H}", "merge", H, None),
    (f"gh pr merge --match-head-commit {H} --merge --repo {P} prv", "merge", H, None),
)  # fmt: skip


@pytest.mark.parametrize(("command", "kind", "sha", "directory"), PRIVATE_FORMS,
                         ids=[c for c, *_ in PRIVATE_FORMS])  # fmt: skip
def test_h22_1_each_private_form_parses_to_a_private_act(
    command: str, kind: str, sha: str | None, directory: str | None
) -> None:
    act = push_grant.granted_act(command)
    assert act is not None, "a private form was not recognised"
    assert (act.kind, act.branch, act.sha, act.directory) == (kind, "prv", sha, directory)
    assert act.remote == "private"


@pytest.mark.parametrize("command", FORMS)
def test_h22_1_each_public_form_is_an_origin_act(command: str) -> None:
    act = push_grant.granted_act(command)
    assert act is not None
    assert act.remote == "origin"


PRIVATE_NOT_FORMS = {
    # §4.2 item 1: the push.
    "push-to-another-remote": f"git push upstream {H}:refs/heads/prv",
    "push-force": f"git push --force private {H}:refs/heads/prv",
    "push-f": f"git push -f private {H}:refs/heads/prv",
    "push-plus-sha": f"git push private +{H}:refs/heads/prv",
    "push-short-sha": f"git push private {SHORT}:refs/heads/prv",
    "push-bare-colon-branch": "git push private :prv",
    "push-bare-colon-ref": "git push private :refs/heads/prv",
    "push-two-refs": f"git push private {H}:refs/heads/prv {H}:refs/heads/other",
    "push-master": f"git push private {H}:refs/heads/master",
    # §4.2 item 2: the repository named.
    "create-other-repo": "gh pr create -R someone/else --head prv",
    "create-repo-twice": f"gh pr create -R {P} -R {P} --head prv",
    "create-repo-no-value": "gh pr create --head prv -R",
    "merge-other-repo": f"gh pr merge prv -R someone/else --merge --match-head-commit {H}",
    "merge-repo-twice": f"gh pr merge prv -R {P} -R {P} --merge --match-head-commit {H}",
    # §4.2 item 3: the private merge.
    "merge-without-merge": f"gh pr merge prv -R {P} --match-head-commit {H}",
    "merge-squash": f"gh pr merge prv -R {P} --merge --squash --match-head-commit {H}",
    "merge-rebase": f"gh pr merge prv -R {P} --merge --rebase --match-head-commit {H}",
    "merge-auto": f"gh pr merge prv -R {P} --merge --auto --match-head-commit {H}",
    "merge-admin": f"gh pr merge prv -R {P} --merge --admin --match-head-commit {H}",
    "merge-d": f"gh pr merge prv -R {P} --merge -d --match-head-commit {H}",
    "merge-delete-branch":
        f"gh pr merge prv -R {P} --merge --delete-branch --match-head-commit {H}",
    "merge-second-positional": f"gh pr merge prv extra -R {P} --merge --match-head-commit {H}",
    "merge-pr-number": f"gh pr merge 12 -R {P} --merge --match-head-commit {H}",
    "merge-without-match-head-commit": f"gh pr merge prv -R {P} --merge",
}  # fmt: skip


@pytest.mark.parametrize("command", PRIVATE_NOT_FORMS.values(), ids=PRIVATE_NOT_FORMS.keys())
def test_h22_2_a_private_command_outside_the_forms_is_not_a_granted_act(command: str) -> None:
    assert push_grant.granted_act(command) is None


@pytest.mark.parametrize(
    "command",
    [f"gh pr merge pub -R expertanalytics/rasputin --match-head-commit {H}",
     "gh pr create -R expertanalytics/rasputin --head pub"],
    ids=["merge", "create"],
)  # fmt: skip
def test_h22_3_the_public_forms_take_no_repository_option(command: str) -> None:
    assert push_grant.granted_act(command) is None


# ---------------------------------------------------------------- Runner

CONFIG_EXITED_128 = (
    "git config could not be read (it exited 128), so it cannot be told whether the "
    "private remote is set up; fix the file git names"
)
CONFIG_TIMED_OUT = (
    "git config did not answer within 10 s, so it cannot be told whether the private "
    "remote is set up"
)
MISSING_REF = (
    "the private remote is set up but refs/remotes/private/master is missing, so it cannot "
    "be told whether this branch holds private commits; run git fetch private"
)
NO_PRIVATE_ONLY = (
    "refs/remotes/private/master holds no commit that origin/master lacks, so it cannot be "
    "told whether this branch holds private commits; run git fetch private"
)


def leak(head: str) -> str:
    return (
        f"head {head[:7]} holds commits of private/master, so it goes only to the "
        "private repository"
    )


def public_only(head: str) -> str:
    return f"head {head[:7]} holds no commit of private/master, so it goes only to origin"


def view(head: str, **change: Any) -> dict[str, Any]:
    """A green `gh pr view --json state,baseRefName,headRefOid,statusCheckRollup`."""
    green = {"state": "OPEN", "baseRefName": "master", "headRefOid": head,
             "statusCheckRollup": [check_run("CI result")]}  # fmt: skip
    return {**green, **change}


def check_run(name: str, status: str = "COMPLETED", conclusion: str = "SUCCESS") -> dict[str, Any]:
    return {"__typename": "CheckRun", "name": name, "status": status, "conclusion": conclusion}


def status_context(name: str, state: str) -> dict[str, Any]:
    return {"__typename": "StatusContext", "context": name, "state": state}


Answer = Callable[[list[str], dict[str, Any]], "subprocess.CompletedProcess[Any]"]


def completed(argv: list[str], kw: dict[str, Any], code: int, out: str) -> Any:
    text = kw.get("text") or kw.get("universal_newlines") or kw.get("encoding")
    return subprocess.CompletedProcess(argv, code, out if text else out.encode(),
                                       "" if text else b"")  # fmt: skip


class Runner:
    """Stands in for `subprocess.run` while a test runs: records every argv, answers
    every `gh` call itself, and the `git config --get remote.private.url` call
    too when `config` is set. Everything else runs for real."""

    def __init__(self, real: Callable[..., Any]) -> None:
        self.real = real
        self.calls: list[tuple[list[str], dict[str, Any]]] = []
        self.gh: Answer = lambda argv, kw: completed(argv, kw, 97, "")
        self.config: Answer | None = None

    def __call__(self, args: Any, *rest: Any, **kw: Any) -> Any:
        argv = [str(word) for word in args]
        self.calls.append((argv, kw))
        if argv[:1] == ["gh"]:
            return self.gh(argv, kw)
        if self.config is not None and "remote.private.url" in argv and "config" in argv:
            return self.config(argv, kw)
        return self.real(args, *rest, **kw)

    def answer_gh(self, found: Any) -> None:
        self.gh = lambda argv, kw: completed(argv, kw, 0, json.dumps(found))

    def ran(self, word: str) -> list[tuple[list[str], dict[str, Any]]]:
        return [(argv, kw) for argv, kw in self.calls if word in argv]


@pytest.fixture
def private(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> PrivateGranted:
    """Fixture A, with `harness_mode.ROOT` pointed at it, so `judge` reads its state,
    never the real checkout's, and a `gh` on PATH that fails if ever reached."""
    built = PrivateGranted(tmp_path.resolve() / "repo")
    monkeypatch.setattr(push_grant.hm, "ROOT", built.repo)
    state = push_grant.hm.state_dir(built.repo)
    assert state is not None and state.resolve().is_relative_to(tmp_path.resolve()), state
    stub = tmp_path / "bin" / "gh"
    stub.parent.mkdir()
    stub.write_text("#!/bin/sh\necho 'the real gh was reached' >&2\nexit 97\n")
    stub.chmod(0o755)
    monkeypatch.setenv("PATH", f"{stub.parent}{os.pathsep}{os.environ['PATH']}")
    for name in ("GIT_DIR", "GIT_WORK_TREE", "GIT_COMMON_DIR", "GIT_INDEX_FILE"):
        monkeypatch.delenv(name, raising=False)
    return built


@pytest.fixture
def runner(private: PrivateGranted, monkeypatch: pytest.MonkeyPatch) -> Runner:
    """Installed after Fixture A is built, so the fixture's own git runs as before."""
    made = Runner(subprocess.run)
    monkeypatch.setattr(subprocess, "run", made)
    return made


def broken_global_config(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """`GIT_CONFIG_GLOBAL` naming a file git cannot parse (exit 128, §4.3)."""
    bad = tmp_path / "broken.gitconfig"
    bad.write_text("[broken\n")
    monkeypatch.setenv("GIT_CONFIG_GLOBAL", str(bad))


# ---------------------------------------------------------------- test 4: home


def test_h22_4a_prv_is_private_and_pub_is_origin(private: PrivateGranted) -> None:
    assert push_grant.home(private.repo, private.prv) == ("private", "")
    assert push_grant.home(private.repo, private.pub) == ("origin", "")


def test_h22_4b_url_set_and_ref_deleted_is_the_missing_ref_reason(private: PrivateGranted) -> None:
    private.delete_private_ref()
    for head in (private.prv, private.pub):
        assert push_grant.home(private.repo, head) == ("", MISSING_REF)


def test_h22_4c_url_set_and_no_private_only_commit_is_the_n0_reason(
    private: PrivateGranted,
) -> None:
    private.private_ref_at_m0()
    for head in (private.prv, private.pub):
        assert push_grant.home(private.repo, head) == ("", NO_PRIVATE_ONLY)


def test_h22_4d_url_unset_and_ref_deleted_is_origin(private: PrivateGranted) -> None:
    private.unset_url()
    private.delete_private_ref()
    for head in (private.prv, private.pub):
        assert push_grant.home(private.repo, head) == ("origin", "")


def test_h22_4d_url_unset_and_ref_left_behind_is_still_read(private: PrivateGranted) -> None:
    """§4.3, the sixth row: an URL removed, a ref left behind: the ref is still read."""
    private.unset_url()
    assert push_grant.home(private.repo, private.prv) == ("private", "")
    assert push_grant.home(private.repo, private.pub) == ("origin", "")


def test_h22_4e_an_unreadable_git_config_is_the_config_reason_and_no_count_runs(
    private: PrivateGranted, runner: Runner, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    broken_global_config(tmp_path, monkeypatch)
    for head in (private.prv, private.pub):
        assert push_grant.home(private.repo, head) == ("", CONFIG_EXITED_128)
    assert runner.ran("remote.private.url"), "home never read remote.private.url"
    assert runner.ran("rev-list") == [], "the counts ran after git config failed"


def test_h22_4f_a_git_config_that_times_out_is_the_timeout_reason(
    private: PrivateGranted, runner: Runner
) -> None:
    """Design review round 3: the config call raises `TimeoutExpired`."""

    def hang(argv: list[str], kw: dict[str, Any]) -> Any:
        raise subprocess.TimeoutExpired(argv, kw.get("timeout") or 10)

    runner.config = hang
    for head in (private.prv, private.pub):
        assert push_grant.home(private.repo, head) == ("", CONFIG_TIMED_OUT)
    [(_, kw), *_] = runner.ran("remote.private.url")
    assert kw.get("timeout") == 10, f"the config call's timeout is {kw.get('timeout')!r}, not 10 s"
    assert runner.ran("rev-list") == [], "the counts ran after git config timed out"


# ---------------------------------------------------------------- test 5: judge


def judged(private: PrivateGranted, command: str) -> Any:
    """`push_grant.judge` on a Bash event in Fixture A, with the reason guard_push gives."""
    reason = push_grant.PR if command.startswith("gh ") else push_grant.PUSH
    return push_grant.judge({"cwd": str(private.repo)}, command, [reason])


def refused_with(grant: Any, reason: str) -> None:
    assert grant is not None, "not judged as a granted act"
    assert not grant.passed, "passed; expected a refusal"
    assert reason in grant.why, grant.why


def used_lines(private: PrivateGranted) -> list[dict[str, Any]]:
    return jsonl(grants_used_path(private.repo))


def test_h22_5a_a_private_push_of_prv_passes_and_logs_the_remote(
    private: PrivateGranted, runner: Runner
) -> None:
    grant = judged(private, f"git push private {private.prv}:refs/heads/prv")
    assert grant is not None and grant.passed, grant
    [line] = used_lines(private)
    assert (line["branch"], line["act"], line["head"], line["remote"]) == (
        "prv", "push", private.prv, "private")  # fmt: skip


def test_h22_5a_an_origin_push_of_pub_passes_and_logs_origin(
    private: PrivateGranted, runner: Runner
) -> None:
    grant = judged(private, f"git push origin {private.pub}:refs/heads/pub")
    assert grant is not None and grant.passed, grant
    [line] = used_lines(private)
    assert line["remote"] == "origin"


def test_h22_5b_an_origin_push_of_prv_is_the_leak_refusal(
    private: PrivateGranted, runner: Runner
) -> None:
    refused_with(judged(private, f"git push origin {private.prv}:refs/heads/prv"),
                 leak(private.prv))  # fmt: skip
    assert used_lines(private) == []


@pytest.mark.parametrize("form", ["create", "merge"])
def test_h22_5b_the_other_origin_forms_of_prv_are_refused_too(
    private: PrivateGranted, runner: Runner, form: str
) -> None:
    command = {"create": "gh pr create --head prv",
               "merge": f"gh pr merge prv --match-head-commit {private.prv}"}[form]  # fmt: skip
    refused_with(judged(private, command), leak(private.prv))
    assert runner.ran("gh") == []


def test_h22_5c_a_private_push_of_pub_goes_only_to_origin(
    private: PrivateGranted, runner: Runner
) -> None:
    refused_with(judged(private, f"git push private {private.pub}:refs/heads/pub"),
                 public_only(private.pub))  # fmt: skip
    assert used_lines(private) == []


def test_h22_5d_origin_push_of_prv_with_the_ref_deleted_is_refused(
    private: PrivateGranted, runner: Runner
) -> None:
    private.delete_private_ref()
    refused_with(judged(private, f"git push origin {private.prv}:refs/heads/prv"), MISSING_REF)


def test_h22_5e_origin_push_of_prv_with_no_private_only_commit_is_refused(
    private: PrivateGranted, runner: Runner
) -> None:
    private.private_ref_at_m0()
    refused_with(judged(private, f"git push origin {private.prv}:refs/heads/prv"),
                 NO_PRIVATE_ONLY)  # fmt: skip


def test_h22_5f_a_private_push_without_the_private_remote_is_refused(
    private: PrivateGranted, runner: Runner
) -> None:
    private.unset_url()
    private.delete_private_ref()
    refused_with(judged(private, f"git push private {private.prv}:refs/heads/prv"),
                 public_only(private.prv))  # fmt: skip


def test_h22_5g_origin_push_of_prv_with_git_config_unreadable_is_refused(
    private: PrivateGranted, runner: Runner
) -> None:
    """The config read alone exits 128: a broken `GIT_CONFIG_GLOBAL` stops every git
    call `judge` makes before `home` (`git rev-parse` exits 128 too), so only the
    one call is made to fail here."""
    runner.config = lambda argv, kw: completed(argv, kw, 128, "")
    refused_with(judged(private, f"git push origin {private.prv}:refs/heads/prv"),
                 CONFIG_EXITED_128)  # fmt: skip
    assert used_lines(private) == []


def test_h22_5g_origin_push_of_prv_with_git_config_timing_out_is_refused(
    private: PrivateGranted, runner: Runner
) -> None:
    def hang(argv: list[str], kw: dict[str, Any]) -> Any:
        raise subprocess.TimeoutExpired(argv, kw.get("timeout") or 10)

    runner.config = hang
    refused_with(judged(private, f"git push origin {private.prv}:refs/heads/prv"),
                 CONFIG_TIMED_OUT)  # fmt: skip


def test_h22_5g_a_broken_global_git_config_never_passes_an_origin_push_of_prv(
    private: PrivateGranted, runner: Runner, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    broken_global_config(tmp_path, monkeypatch)
    try:
        grant = judged(private, f"git push origin {private.prv}:refs/heads/prv")
    except subprocess.CalledProcessError:
        return  # the guard then fails closed: an exception in the hook is no pass
    assert grant is None or not grant.passed, grant
    assert used_lines(private) == []


# ---------------------------------------------------------------- test 6: the walk


def advance_private_master(private: PrivateGranted, relative: str) -> str:
    """A code commit on private/master (P2), made on a scratch branch then deleted."""
    git(private.repo, "checkout", "-q", "-b", "scratch-private", "refs/remotes/private/master")
    sha = commit_file(private.repo, relative, "p = 2\n", f"P2 changes {relative}")
    git(private.repo, "update-ref", "refs/remotes/private/master", sha)
    git(private.repo, "checkout", "-q", "prv")
    git(private.repo, "branch", "-q", "-D", "scratch-private")
    return sha


def merge_into_prv(private: PrivateGranted, ref: str) -> str:
    git(private.repo, "merge", "-q", "--no-ff", "--no-edit", ref)
    return private.head("prv")


def walk(private: PrivateGranted) -> str | None:
    found: str | None = push_grant.reviewed(
        private.repo, private.head("prv"), private.records(), remote="private"
    )
    return found


def test_h22_6_a_code_commit_on_prv_after_the_review_is_refused(private: PrivateGranted) -> None:
    git(private.repo, "checkout", "-q", "prv")
    late = commit_file(private.repo, "tools/y.py", "y = 1\n", "a late code change")
    refuses_naming(walk(private), late, "tools/y.py", private.prv)


def test_h22_6_a_clean_merge_of_private_master_after_the_review_passes(
    private: PrivateGranted,
) -> None:
    git(private.repo, "checkout", "-q", "prv")
    advance_private_master(private, "src/p2.py")
    merge_into_prv(private, "refs/remotes/private/master")
    assert walk(private) is None


def test_h22_6_a_clean_merge_of_origin_master_into_prv_is_refused(
    private: PrivateGranted,
) -> None:
    git(private.repo, "checkout", "-q", "prv")
    merge = merge_into_prv(private, "refs/remotes/origin/master")
    found = walk(private)
    assert found is not None, "a merge of origin/master passed as a merge of private/master"
    assert f"commit {merge[:7]}" in found, found


# ---------------------------------------------------------------- test 7: checks

HEAD2 = "fedcba9876543210fedcba9876543210fedcba98"
UNREAD = "the pull request's checks could not be read"

CHECKS = {
    "state-merged": ({"state": "MERGED"}, "the pull request is MERGED, not open"),
    "base-main": ({"baseRefName": "main"}, "the pull request is into main, not master"),
    "head-other": ({"headRefOid": HEAD2}, f"the pull request's head is {HEAD2[:7]}, not {H[:7]}"),
    "rollup-empty": ({"statusCheckRollup": []}, "no checks reported"),
    "checkrun-in-progress": (
        {"statusCheckRollup": [check_run("ok"), check_run("build", "IN_PROGRESS", "")]}, "build"),
    "checkrun-failure": (
        {"statusCheckRollup": [check_run("ok"), check_run("lint", conclusion="FAILURE")]}, "lint"),
    "statuscontext-pending": (
        {"statusCheckRollup": [check_run("ok"), status_context("ci/x", "PENDING")]}, "ci/x"),
    "all-skipped": (
        {"statusCheckRollup": [check_run("a", conclusion="SKIPPED"),
                               check_run("b", conclusion="NEUTRAL")]}, ""),
}  # fmt: skip


def test_h22_7_a_green_view_passes() -> None:
    assert push_grant.checks(view(H), H) is None


def test_h22_7_skipped_and_neutral_beside_a_success_pass() -> None:
    rollup = [check_run("a"), check_run("b", conclusion="SKIPPED"),
              check_run("c", conclusion="NEUTRAL"), status_context("d", "SUCCESS")]  # fmt: skip
    assert push_grant.checks(view(H, statusCheckRollup=rollup), H) is None


@pytest.mark.parametrize(("change", "reason"), CHECKS.values(), ids=CHECKS.keys())
def test_h22_7_each_condition_refuses_with_its_own_reason(
    change: dict[str, Any], reason: str
) -> None:
    found = push_grant.checks(view(H, **change), H)
    assert found is not None, "passed; expected a refusal"
    assert reason in found, found


@pytest.mark.parametrize("found", [[], {}, "text", None], ids=["list", "empty", "str", "null"])
def test_h22_7_json_without_the_keys_cannot_be_read(found: Any) -> None:
    assert push_grant.checks(found, H) == UNREAD


MALFORMED = {
    "item-a-string": ["CI result"],
    "item-null": [None],
    "item-a-list-beside-a-green-run": [{"__typename": "CheckRun", "name": "ok",
                                        "status": "COMPLETED", "conclusion": "SUCCESS"}, []],
    "red-checkrun-without-a-name": [{"__typename": "CheckRun", "status": "COMPLETED",
                                     "conclusion": "FAILURE"}],
    "red-checkrun-without-a-state": [{"__typename": "CheckRun", "name": "lint"}],
    "red-statuscontext-without-a-context": [{"__typename": "StatusContext", "state": "FAILURE"}],
    "red-statuscontext-without-a-state": [{"__typename": "StatusContext", "context": "ci/x"}],
}  # fmt: skip


@pytest.mark.parametrize("rollup", MALFORMED.values(), ids=MALFORMED.keys())
def test_h22_7_a_malformed_check_item_cannot_be_read(rollup: list[Any]) -> None:
    """Review round 1 follow-up: never "check None is None"."""
    assert push_grant.checks(view(H, statusCheckRollup=rollup), H) == UNREAD


def private_merge(private: PrivateGranted) -> str:
    return f"gh pr merge prv -R {P} --merge --match-head-commit {private.prv}"


def test_h22_7_a_private_merge_with_green_checks_passes_and_reads_them_once(
    private: PrivateGranted, runner: Runner
) -> None:
    runner.answer_gh(view(private.prv))
    grant = judged(private, private_merge(private))
    assert grant is not None and grant.passed, grant
    [(argv, kw)] = runner.ran("gh")
    assert argv[:4] == ["gh", "pr", "view", "prv"] and P in argv, argv
    assert kw.get("timeout") == 20, kw
    [line] = used_lines(private)
    assert (line["act"], line["remote"]) == ("merge", "private")


def test_h22_7_a_private_merge_with_red_checks_is_refused(
    private: PrivateGranted, runner: Runner
) -> None:
    runner.answer_gh(view(private.prv, statusCheckRollup=[check_run("lint", conclusion="FAILURE")]))
    refused_with(judged(private, private_merge(private)), "lint")
    assert used_lines(private) == []


def gh_exits_1(argv: list[str], kw: dict[str, Any]) -> Any:
    return completed(argv, kw, 1, "")


def gh_times_out(argv: list[str], kw: dict[str, Any]) -> Any:
    raise subprocess.TimeoutExpired(argv, kw.get("timeout") or 20)


def gh_prints_a_list(argv: list[str], kw: dict[str, Any]) -> Any:
    return completed(argv, kw, 0, "[]")


@pytest.mark.parametrize("answer", [gh_exits_1, gh_times_out, gh_prints_a_list],
                         ids=["exits-1", "times-out", "prints-a-list"])  # fmt: skip
def test_h22_7_checks_that_cannot_be_read_refuse_the_private_merge(
    private: PrivateGranted, runner: Runner, answer: Answer
) -> None:
    runner.gh = answer
    refused_with(judged(private, private_merge(private)), UNREAD)
    assert used_lines(private) == []


# ---------------------------------------------------------------- test 8


def test_h22_8_a_second_private_merge_of_the_same_head_is_refused(
    private: PrivateGranted, runner: Runner
) -> None:
    runner.answer_gh(view(private.prv))
    first = judged(private, private_merge(private))
    assert first is not None and first.passed, first
    second = judged(private, private_merge(private))
    assert second is not None and not second.passed, second
    assert [line["act"] for line in used_lines(private)] == ["merge"]
    assert second.why == f"head {private.prv[:7]} was merged already in this window", second.why


def test_h22_8_a_second_enqueue_of_the_same_head_is_refused(
    private: PrivateGranted, runner: Runner
) -> None:
    command = f"gh pr merge pub --match-head-commit {private.pub}"
    first = judged(private, command)
    assert first is not None and first.passed, first
    second = judged(private, command)
    assert second is not None and not second.passed, second
    assert second.why == f"head {private.pub[:7]} was enqueued already in this window", second.why
    assert [line["act"] for line in used_lines(private)] == ["enqueue"]
    assert runner.ran("gh") == [], "an enqueue read the pull request's checks"


# ================================================================ h24 item 1
#
# docs/increments/h24-harness-fixes.md §3 item 1: `@perf` runs after the code
# review, so its evidence commit always follows the verdict. A post-review
# commit whose only changes beyond prose are files it *adds* under a dated
# `docs/benchmarks/<YYYY-MM-DD>/` folder counts like prose. Each case runs on
# the `built` fixture (review APPROVED at `x`, prose `rec` after it).
# Pinned here beyond the design's words: a deleted dated file is refused like
# an edited one (the design passes only what the commit adds), and the dated
# folder is matched from the start of the path, so `tools/docs/benchmarks/…`
# is code.

DATED = "docs/benchmarks/2026-10-17/x"


def commit_files(repo: Path, files: dict[str, str], message: str) -> str:
    """Write each file (creating its directory), commit them together; the full sha."""
    for relative, text in files.items():
        path = repo / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text)
        git(repo, "add", relative)
    git(repo, "commit", "-q", "-m", message)
    return git(repo, "rev-parse", "HEAD").strip()


def test_h24_1_added_dated_benchmark_evidence_after_the_review_passes(built: Granted) -> None:
    commit_files(built.repo, {f"{DATED}/times.tsv": "n\tt\n1\t2\n",
                              f"{DATED}/summary.json": "{}\n",
                              f"{DATED}/plot.py": "print(1)\n"},
                 "@perf evidence")  # fmt: skip
    assert judge(built) is None


def test_h24_1_an_added_run_json_passes_the_accepted_risk(built: Granted) -> None:
    """A forged run.json could become a later baseline; §3 item 1 accepts that."""
    commit_file(built.repo, f"{DATED}/run.json", '{"wall": 1.0}\n', "@perf run record")
    assert judge(built) is None


def test_h24_1_an_edit_of_an_existing_dated_file_is_refused_naming_it(built: Granted) -> None:
    path = f"{DATED}/summary.json"
    commit_file(built.repo, path, "{}\n", "@perf evidence")
    edit = commit_file(built.repo, path, '{"forged": true}\n', "edit the evidence")
    refuses_naming(judge(built), edit, path, built.x)


def test_h24_1_a_deleted_dated_file_is_refused_naming_it(built: Granted) -> None:
    path = f"{DATED}/summary.json"
    commit_file(built.repo, path, "{}\n", "@perf evidence")
    git(built.repo, "rm", "-q", path)
    git(built.repo, "commit", "-q", "-m", "drop the evidence")
    drop = built.head
    refuses_naming(judge(built), drop, path, built.x)


@pytest.mark.parametrize(
    "path",
    ["docs/benchmarks/quick/new.json", "docs/benchmarks/2026-10-17-x/new.json",
     "docs/benchmarks/new.json", "tools/docs/benchmarks/2026-10-17/new.py"],
    ids=["quick", "not-a-date-folder", "undated", "not-at-the-start"],
)  # fmt: skip
def test_h24_1_an_added_file_outside_a_dated_folder_is_refused(built: Granted, path: str) -> None:
    late = commit_file(built.repo, path, "{}\n", f"add {path}")
    refuses_naming(judge(built), late, path, built.x)


def test_h24_1_dated_evidence_with_a_code_change_is_refused_naming_the_code(
    built: Granted,
) -> None:
    late = commit_files(built.repo, {f"{DATED}/run.json": "{}\n", "tools/x.py": "x = 9\n"},
                        "evidence and a code change")  # fmt: skip
    refuses_naming(judge(built), late, "tools/x.py", built.x)
