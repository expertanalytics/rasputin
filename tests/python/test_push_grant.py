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

import re
from pathlib import Path

import pytest

from harness_fixtures import (
    Granted,
    Tool,
    commit_file,
    git,
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
