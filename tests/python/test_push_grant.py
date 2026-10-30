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
    add_remote,
    commit_file,
    flag_granting,
    flag_path,
    git,
    grants_used_path,
    jsonl,
    make_repo,
    pin_remotes,
    publish,
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


# ---------------------------------------------------------------- redo in form
#
# Ola's ruling of 2026-10-27 (REQUIRED-READING.md, Unattended mode): a granted
# act refused only for its form may be redone in the form its reason names,
# and the refusal says so in the words of guard_push.py's AMEND_FORM.

REDO_IN_FORM = "; it may be redone in that form"


def test_not_bare_says_it_may_be_redone_in_that_form() -> None:
    assert push_grant.NOT_BARE.endswith(f"run it alone{REDO_IN_FORM}"), push_grant.NOT_BARE


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
    add_remote(repo, "origin")
    publish(repo, "origin", "master")
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


# ---------------------------------------------------------------- merge needs a verdict
#
# PR #400 merged head f29480e9, a test fix made after the last approved review.
# A full-lane head merges only when reviewed; after the verdict, only commits
# touching the increment file and ROADMAP.md (docs/increments/README.md step 4).


def merge_why(built: Granted) -> str | None:
    found: str | None = push_grant.merge_refusal(built.repo, built.head, built.records())
    return found


def test_merge_of_the_reviewed_head_and_its_record_passes(built: Granted) -> None:
    commit_file(built.repo, "ROADMAP.md", "row\n", "h99 roadmap row")
    assert merge_why(built) is None


def test_merge_after_a_test_fix_past_the_verdict_is_refused(built: Granted) -> None:
    late = commit_file(built.repo, "tests/python/test_x.py", "t = 1\n", "a test fix")
    refuses_naming(merge_why(built), late, "tests/python/test_x.py", built.x)


def test_merge_after_other_prose_past_the_verdict_is_refused(built: Granted) -> None:
    late = commit_file(built.repo, "docs/notes.md", "n\n", "a prose change")
    refuses_naming(merge_why(built), late, "docs/notes.md", built.x)


def test_merge_without_a_verdict_is_refused(built: Granted) -> None:
    write_reviews(built.repo, [])
    assert NO_REVIEW in (merge_why(built) or "")


def test_merge_of_a_fast_lane_branch_needs_no_verdict(built: Granted) -> None:
    git(built.repo, "checkout", "-q", "-b", "quick", "master")
    commit_file(built.repo, "tools/q.py", "q = 1\n", "a fast-lane change")
    write_reviews(built.repo, [])
    head = git(built.repo, "rev-parse", "HEAD").strip()
    assert push_grant.merge_refusal(built.repo, head, []) is None


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
    pin_remotes(monkeypatch, push_grant, built.repo)
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


def test_h22_5d_origin_push_of_prv_with_the_local_ref_deleted_is_still_the_leak(
    private: PrivateGranted, runner: Runner
) -> None:
    """The private master is read from the private remote, so the local ref's absence
    changes nothing (docs/increments/h38-fastlane-push.md §3, change 4)."""
    private.delete_private_ref()
    refused_with(judged(private, f"git push origin {private.prv}:refs/heads/prv"),
                 leak(private.prv))  # fmt: skip


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
    publish(private.repo, "private", sha)
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


# ---------------------------------------------------------------- h29 §3.2
#
# docs/increments/h29-review-rules.md §4 tests 2 to 4: a grant refused only
# because the head is not reviewed says it may be redone once an approving
# verdict on that head is recorded; any other refusal does not say so. Each
# case judges in process on `Granted` with the flag on and `harness_mode.ROOT`
# pointed at it, as the `private` fixture does, for the three granted forms.


def redo(head: str) -> str:
    return f"may be redone once @reviewer's verdict approving {head[:7]} is recorded"


@pytest.fixture
def grant(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Granted:
    built = Granted(tmp_path.resolve() / "repo")
    monkeypatch.setattr(push_grant.hm, "ROOT", built.repo)
    pin_remotes(monkeypatch, push_grant, built.repo)
    state = push_grant.hm.state_dir(built.repo)
    assert state is not None and state.resolve().is_relative_to(tmp_path.resolve()), state
    for name in ("GIT_DIR", "GIT_WORK_TREE", "GIT_COMMON_DIR", "GIT_INDEX_FILE"):
        monkeypatch.delenv(name, raising=False)
    return built


def acts(head: str) -> dict[str, str]:
    """The three granted forms for `feat` at `head`: push, pull request, enqueue."""
    return {
        "push": f"git push origin {head}:refs/heads/feat",
        "create": "gh pr create --head feat --title 't' --body-file f",
        "enqueue": f"gh pr merge feat --match-head-commit {head}",
    }


def granted_judge(built: Granted, command: str) -> Any:
    reason = push_grant.PR if command.startswith("gh ") else push_grant.PUSH
    return push_grant.judge({"cwd": str(built.repo)}, command, [reason])


@pytest.mark.parametrize("act", ["push", "create", "enqueue"])
def test_h29_2_no_record_is_refused_saying_it_may_be_redone(grant: Granted, act: str) -> None:
    write_reviews(grant.repo, [])
    found = granted_judge(grant, acts(grant.head)[act])
    refused_with(found, NO_REVIEW)
    assert found.why.startswith(f"head {grant.head[:7]} is not reviewed: "), found.why
    assert redo(grant.head) in found.why, found.why


@pytest.mark.parametrize("act", ["push", "create", "enqueue"])
def test_h29_3_changes_requested_on_the_head_is_refused_saying_it_may_be_redone(
    grant: Granted, act: str
) -> None:
    write_reviews(grant.repo, [review_record(grant.head, "CHANGES REQUESTED")])
    found = granted_judge(grant, acts(grant.head)[act])
    refused_with(found, CHANGES)
    assert redo(grant.head) in found.why, found.why


def test_h29_4_an_enqueue_already_used_in_this_window_does_not_say_redo(grant: Granted) -> None:
    command = acts(grant.head)["enqueue"]
    first = granted_judge(grant, command)
    assert first is not None and first.passed, first
    second = granted_judge(grant, command)
    refused_with(second, f"head {grant.head[:7]} was enqueued already in this window")
    assert "redone" not in second.why, second.why


# ---------------------------------------------------------------- the fast lane
#
# docs/increments/h38-fastlane-push.md §3 and §6. Present or away, a branch
# for origin whose diff against origin's real master (read from the remote,
# never from the local tracking ref) changes no full-lane path, no increment
# file, no governed path and no CI-config path, at most 700 net lines, may be
# pushed, opened (when origin's copy of the branch is the local head) and
# enqueued once per head (when its pull request's checks are green), without
# a grant or a review record. Any other branch keeps the grant and review
# rule. Each case judges in process; "present" is the flag file removed.
#
# Pinned here beyond the design's words: `pr_view(where, branch, repo)` names
# its third parameter `repo`; a create's refusal reads "origin has no
# refs/heads/<b>" or "origin's refs/heads/<b> is <sha7>"; a redirected URL's
# refusal names the URL git printed.

PUBLIC_REPO = "expertanalytics/rasputin"
UNREAD = "the pull request's checks could not be read"
MODES = ("present", "away")


class Fast:
    """`master` (1 commit) pushed to a bare `origin`; branch `quick` checked out with
    one commit of `files`, pushed to origin's `quick` too unless `pushed` is false.
    The flag is on and grants nothing. `pr_view` answers `view` and records each
    repository it is asked about in `seen`."""

    def __init__(self, root: Path, files: dict[str, str], pushed: bool = True) -> None:
        self.repo = make_repo(root)
        self.base = git(self.repo, "rev-parse", "master").strip()
        self.origin = add_remote(self.repo, "origin")
        publish(self.repo, "origin", self.base)
        git(self.repo, "checkout", "-q", "-b", "quick")
        self.head = commit_files(self.repo, files, "a fast-lane change")
        if pushed:
            publish(self.repo, "origin", self.head, "quick")
        self.since = flag_granting(self.repo, [])
        self.view: Any = view(self.head)
        self.seen: list[str] = []

    def pr_view(self, where: str, branch: str, repo: str = "") -> Any:
        self.seen.append(repo)
        return self.view

    def acts(self) -> dict[str, str]:
        return {
            "push": f"git push origin {self.head}:refs/heads/quick",
            "create": "gh pr create --head quick --title 't' --body-file f",
            "enqueue": f"gh pr merge quick --match-head-commit {self.head}",
        }

    def judge(self, act: str) -> Any:
        command = self.acts()[act]
        reason = push_grant.PR if command.startswith("gh ") else push_grant.PUSH
        return push_grant.judge({"cwd": str(self.repo)}, command, [reason])

    def in_mode(self, mode: str) -> None:
        if mode == "present":
            flag_path(self.repo).unlink()

    def used(self) -> list[Any]:
        return jsonl(grants_used_path(self.repo))


def lines_of(n: int) -> str:
    return "".join(f"v{i} = {i}\n" for i in range(n))


FastRepo = Callable[..., Fast]


@pytest.fixture
def fast_repo(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> FastRepo:
    for name in ("GIT_DIR", "GIT_WORK_TREE", "GIT_COMMON_DIR", "GIT_INDEX_FILE", "GH_REPO"):
        monkeypatch.delenv(name, raising=False)

    def build(files: dict[str, str], pushed: bool = True) -> Fast:
        made = Fast(tmp_path.resolve() / "repo", files, pushed)
        monkeypatch.setattr(push_grant.hm, "ROOT", made.repo)
        pin_remotes(monkeypatch, push_grant, made.repo)
        monkeypatch.setitem(push_grant.judge.__globals__, "pr_view", made.pr_view)
        return made

    return build


def not_judged_or_refused(found: Any, *words: str) -> None:
    """None, or a refusal whose reason holds each of `words`."""
    if found is not None:
        assert not found.passed, f"passed: {found}"
        for word in words:
            assert word in found.why, found.why


@pytest.mark.parametrize("act", ["push", "create", "enqueue"])
def test_a_fast_lane_branch_passes_without_grant_or_review(fast_repo: FastRepo, act: str) -> None:
    fast = fast_repo({"tools/y.py": "y = 1\n"})
    found = fast.judge(act)
    assert found is not None and found.passed, found
    [line] = fast.used()
    assert (line["branch"], line["head"], line["lane"]) == ("quick", fast.head, "fast"), line


def test_a_fast_lane_branch_of_exactly_700_lines_passes(fast_repo: FastRepo) -> None:
    found = fast_repo({"tools/big.py": lines_of(700)}).judge("push")
    assert found is not None and found.passed, found


def test_a_fast_lane_branch_is_enqueued_once_per_head(fast_repo: FastRepo) -> None:
    fast = fast_repo({"tools/y.py": "y = 1\n"})
    assert fast.judge("enqueue").passed
    refused_with(fast.judge("enqueue"), "enqueued already in this window")


# §6 test 1 -- present: push and create pass, logged with no window.


@pytest.mark.parametrize("act", ["push", "create"])
def test_present_a_fast_lane_push_and_create_pass_and_log_no_window(
    fast_repo: FastRepo, act: str
) -> None:
    fast = fast_repo({"tools/y.py": "y = 1\n"})
    fast.in_mode("present")
    found = fast.judge(act)
    assert found is not None and found.passed, found
    [line] = fast.used()
    assert (line["lane"], line["window"], line["head"]) == ("fast", None, fast.head), line


# §6 test 2 -- present: a green enqueue passes, read from the public repository, once.


def test_present_a_green_fast_lane_enqueue_passes_once_per_head(fast_repo: FastRepo) -> None:
    fast = fast_repo({"tools/y.py": "y = 1\n"})
    fast.in_mode("present")
    found = fast.judge("enqueue")
    assert found is not None and found.passed, found
    assert fast.seen and set(fast.seen) == {PUBLIC_REPO}, fast.seen
    refused_with(fast.judge("enqueue"), "enqueued already")
    assert len(fast.used()) == 1


# §6 test 3 -- CI not green: the enqueue is refused, the reason names why.

NOT_GREEN: dict[str, tuple[Callable[[str], Any], Callable[[str], str]]] = {
    "pending": (lambda h: view(h, statusCheckRollup=[check_run("CI result", "IN_PROGRESS", "")]),
                lambda h: "check CI result is IN_PROGRESS"),
    "failed": (lambda h: view(h, statusCheckRollup=[check_run("CI result", conclusion="FAILURE")]),
               lambda h: "check CI result is FAILURE"),
    "another-head": (lambda h: view(H),
                     lambda h: f"the pull request's head is {H[:7]}, not {h[:7]}"),
    "none-reported": (lambda h: view(h, statusCheckRollup=[]), lambda h: "no checks reported"),
    "unreadable": (lambda h: None, lambda h: UNREAD),
}  # fmt: skip


@pytest.mark.parametrize("case", NOT_GREEN)
@pytest.mark.parametrize("mode", MODES)
def test_a_fast_lane_enqueue_without_green_checks_is_refused(
    fast_repo: FastRepo, mode: str, case: str
) -> None:
    fast = fast_repo({"tools/y.py": "y = 1\n"})
    fast.in_mode(mode)
    answer, reason = NOT_GREEN[case]
    fast.view = answer(fast.head)
    refused_with(fast.judge("enqueue"), reason(fast.head))
    assert set(fast.seen) == {PUBLIC_REPO}, fast.seen
    assert fast.used() == []


# §6 test 4 -- a governed or CI-config path keeps the branch out of the fast lane.

GOVERNED_OR_CI = (
    "CLAUDE.md", ".claude/settings.json", ".claude/agents/x.md", ".claude/hooks/x.py",
    "tools/check_x.py", "tools/push_grant.py", ".github/workflows/main.yaml", "pyproject.toml",
    "tests/python/conftest.py", "tests/cpp/CMakeLists.txt", "pytest.ini", "ruff.toml",
)  # fmt: skip


@pytest.mark.parametrize("act", ["push", "create", "enqueue"])
@pytest.mark.parametrize("path", GOVERNED_OR_CI)
@pytest.mark.parametrize("mode", MODES)
def test_a_branch_changing_a_governed_or_ci_config_path_is_not_judged(
    fast_repo: FastRepo, mode: str, path: str, act: str
) -> None:
    fast = fast_repo({"tools/y.py": "y = 1\n", path: "changed = 1\n"})
    fast.in_mode(mode)
    assert fast.judge(act) is None
    assert fast.used() == []


# §6 test 5 -- the other cases outside the lane, in both modes.


@pytest.mark.parametrize(
    "files",
    [{"tools/y.py": "y = 1\n", "src/a.cpp": "int a;\n"},
     {"include/terrain/mesh/m.hpp": "int m;\n"},
     {"src_python/tin_engine/crs.py": "c = 1\n"},
     {"tools/y.py": "y = 1\n", "docs/increments/h99-x.md": "# x\n"},
     {"tools/big.py": lines_of(701)}],
    ids=["full-lane-src", "full-lane-mesh", "full-lane-crs", "increment-file", "over-700"],
)  # fmt: skip
@pytest.mark.parametrize("act", ["push", "create", "enqueue"])
@pytest.mark.parametrize("mode", MODES)
def test_any_other_ungranted_branch_is_not_judged(
    fast_repo: FastRepo, mode: str, files: dict[str, str], act: str
) -> None:
    fast = fast_repo(files)
    fast.in_mode(mode)
    assert fast.judge(act) is None
    assert fast.used() == []


@pytest.mark.parametrize("mode", MODES)
def test_a_fast_lane_push_naming_another_sha_is_not_judged(fast_repo: FastRepo, mode: str) -> None:
    fast = fast_repo({"tools/y.py": "y = 1\n"})
    fast.in_mode(mode)
    command = f"git push origin {H}:refs/heads/quick"
    assert push_grant.judge({"cwd": str(fast.repo)}, command, [push_grant.PUSH]) is None


# §6 test 6 -- the remotes' real master, never the local tracking refs.


def forge(repo: Path, ref: str, sha: str) -> None:
    """Write `ref` through its file, as a plain shell write would, past every git command."""
    (repo / ".git" / ref).write_text(sha + "\n")
    assert git(repo, "rev-parse", ref).strip() == sha


@pytest.mark.parametrize("mode", MODES)
def test_a_forged_local_origin_master_does_not_open_the_fast_lane(
    fast_repo: FastRepo, mode: str
) -> None:
    fast = fast_repo({"src/x.cpp": "int x;\n"})
    fast.in_mode(mode)
    forge(fast.repo, "refs/remotes/origin/master", fast.head)
    for act in ("push", "create", "enqueue"):
        assert fast.judge(act) is None, act
    assert fast.used() == []


def test_a_forged_local_private_master_does_not_send_private_commits_to_origin(
    private: PrivateGranted, runner: Runner
) -> None:
    tree = git(private.repo, "rev-parse", f"{private.m0}^{{tree}}").strip()
    unrelated = git(private.repo, "commit-tree", tree, "-m", "an unrelated history").strip()
    forge(private.repo, "refs/remotes/private/master", unrelated)
    refused_with(judged(private, f"git push origin {private.prv}:refs/heads/prv"),
                 "holds commits of private/master")  # fmt: skip
    assert used_lines(private) == []


@pytest.mark.parametrize("mode", MODES)
def test_a_stale_local_origin_master_is_not_what_the_lane_is_judged_against(
    fast_repo: FastRepo, mode: str
) -> None:
    """origin's master moved to the branch's full-lane commit, the tracking ref not
    told: against the real master the branch changes only `tools/y.py`."""
    fast = fast_repo({"src/a.cpp": "int a;\n"})
    full = fast.head
    fast.head = commit_files(fast.repo, {"tools/y.py": "y = 1\n"}, "a fast-lane change on top")
    publish(fast.repo, "origin", fast.head, "quick")
    git(fast.origin, "update-ref", "refs/heads/master", full)  # `quick` brought its objects
    assert git(fast.repo, "rev-parse", "refs/remotes/origin/master").strip() == fast.base
    fast.in_mode(mode)
    found = fast.judge("push")
    assert found is not None and found.passed, found


@pytest.mark.parametrize("mode", MODES)
def test_an_unfetched_origin_master_is_fetched_and_judged_against(
    fast_repo: FastRepo, mode: str, tmp_path: Path
) -> None:
    fast = fast_repo({"tools/y.py": "y = 1\n"})
    other = tmp_path / "other"
    git(tmp_path, "clone", "-q", str(fast.origin), str(other))
    git(other, "config", "user.email", "t@example.invalid")
    git(other, "config", "user.name", "T")
    moved = commit_file(other, "src/b.cpp", "int b;\n", "master moves on origin")
    git(other, "push", "-q", "origin", "HEAD:refs/heads/master")
    missing = subprocess.run(["git", "-C", str(fast.repo), "cat-file", "-e", moved],
                             capture_output=True, check=False)  # fmt: skip
    assert missing.returncode != 0, "the moved master is already here"
    fast.in_mode(mode)
    found = fast.judge("push")
    assert found is not None and found.passed, found
    assert git(fast.repo, "cat-file", "-t", moved).strip() == "commit"


# §6 test 7 -- a create needs origin's copy of the branch at the local head.


@pytest.mark.parametrize("mode", MODES)
def test_a_create_with_no_branch_on_origin_is_refused(fast_repo: FastRepo, mode: str) -> None:
    fast = fast_repo({"tools/y.py": "y = 1\n"}, pushed=False)
    fast.in_mode(mode)
    refused_with(fast.judge("create"), "origin has no refs/heads/quick")
    assert fast.used() == []


@pytest.mark.parametrize("mode", MODES)
def test_a_create_with_origins_branch_elsewhere_is_refused(fast_repo: FastRepo, mode: str) -> None:
    fast = fast_repo({"tools/y.py": "y = 1\n"})
    publish(fast.repo, "origin", fast.base, "quick")
    fast.in_mode(mode)
    refused_with(fast.judge("create"), f"origin's refs/heads/quick is {fast.base[:7]}")
    assert fast.used() == []


# §6 tests 8 and 9 -- an unreadable remote, a broken flag.


@pytest.mark.parametrize("mode", MODES)
def test_an_origin_that_cannot_be_read_passes_nothing(
    fast_repo: FastRepo, mode: str, tmp_path: Path
) -> None:
    fast = fast_repo({"tools/y.py": "y = 1\n"})
    gone = str(tmp_path / "gone.git")
    git(fast.repo, "remote", "set-url", "origin", gone)
    push_grant.judge.__globals__["REMOTE_URLS"]["origin"] = gone
    fast.in_mode(mode)
    for act in ("push", "create", "enqueue"):
        not_judged_or_refused(fast.judge(act))
    assert fast.used() == []


def test_with_the_flag_broken_a_fast_lane_branch_is_not_judged(fast_repo: FastRepo) -> None:
    fast = fast_repo({"tools/y.py": "y = 1\n"})
    flag_path(fast.repo).write_text("{not json")
    for act in ("push", "create", "enqueue"):
        assert fast.judge(act) is None, act
    assert fast.used() == []


# §6 test 11 -- a fetch or push URL that is not the pinned one, or gh pointed elsewhere.


def second_origin(fast: Fast, tmp_path: Path) -> str:
    """A second bare repository holding origin's master: a redirect's target."""
    second = tmp_path / "second.git"
    git(tmp_path, "init", "-q", "--bare", str(second))
    git(fast.repo, "push", "-q", str(second), f"{fast.base}:refs/heads/master")
    return str(second)


REDIRECTS = ("pushInsteadOf", "insteadOf", "GIT_CONFIG_COUNT", "added-pushurl", "second-url",
             "gh-resolved", "GH_REPO")  # fmt: skip


@pytest.mark.parametrize("act", ["push", "create", "enqueue"])
@pytest.mark.parametrize("redirect", REDIRECTS)
@pytest.mark.parametrize("mode", MODES)
def test_a_redirected_remote_passes_nothing(
    fast_repo: FastRepo, tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
    mode: str, redirect: str, act: str,
) -> None:  # fmt: skip
    fast = fast_repo({"tools/y.py": "y = 1\n"})
    second, origin = second_origin(fast, tmp_path), str(fast.origin)
    named = {"gh-resolved": "gh-resolved", "GH_REPO": "GH_REPO"}.get(redirect, second)
    if redirect in ("pushInsteadOf", "insteadOf"):
        config = tmp_path / "global.gitconfig"
        config.write_text(f'[url "{second}"]\n\t{redirect} = {origin}\n')
        monkeypatch.setenv("GIT_CONFIG_GLOBAL", str(config))
    elif redirect == "GIT_CONFIG_COUNT":
        monkeypatch.setenv("GIT_CONFIG_COUNT", "1")
        monkeypatch.setenv("GIT_CONFIG_KEY_0", f"url.{second}.insteadOf")
        monkeypatch.setenv("GIT_CONFIG_VALUE_0", origin)
    elif redirect == "added-pushurl":
        git(fast.repo, "config", "remote.origin.pushurl", second)
    elif redirect == "second-url":
        git(fast.repo, "config", "--add", "remote.origin.url", second)
    elif redirect == "gh-resolved":
        git(fast.repo, "config", "remote.private.gh-resolved", "base")
    else:
        monkeypatch.setenv("GH_REPO", "someone/else")
    fast.in_mode(mode)
    not_judged_or_refused(fast.judge(act), named)
    assert fast.used() == []


def test_a_redirected_private_push_url_refuses_a_granted_push_to_origin(
    private: PrivateGranted, runner: Runner, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    elsewhere = tmp_path / "elsewhere.git"
    git(tmp_path, "init", "-q", "--bare", str(elsewhere))
    config = tmp_path / "global.gitconfig"
    config.write_text(f'[url "{elsewhere}"]\n\tpushInsteadOf = {private.url}\n')
    monkeypatch.setenv("GIT_CONFIG_GLOBAL", str(config))
    found = judged(private, f"git push origin {private.pub}:refs/heads/pub")
    refused_with(found, str(elsewhere))
    assert "private" in found.why, found.why
    assert used_lines(private) == []


def test_a_granted_full_lane_branch_still_needs_its_review(grant: Granted) -> None:
    write_reviews(grant.repo, [])
    refused_with(granted_judge(grant, acts(grant.head)["push"]), NO_REVIEW)
