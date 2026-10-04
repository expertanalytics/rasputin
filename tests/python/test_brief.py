"""`tools/brief.py`: the fixed part of a brief (h9, tests 1 to 9).

Spec: `docs/increments/h9-spawn-briefs.md` §3.1, §3.1a, §3.2 and §3.3. The
script is run from the fixture's copy (`brief_fixtures.make_brief_repo`), from the
main checkout, with a temporary `HOME` holding the transcript and
`CLAUDE_CODE_SESSION_ID` set; linked worktrees A, B and C stand beside it. The
block's format is checked against `brief_fixtures`, which writes §3.2 a second
time from the design's text.

One in-process test (read-only is derived) loads the fixture's copy and calls
its `main(argv)`, so `WRITES` can be changed for the run; see the handback.
"""

from __future__ import annotations

import re
from datetime import datetime
from pathlib import Path
from typing import Any

import pytest

from brief_fixtures import (
    BRIEF,
    COMMON,
    NEEDS_INCREMENT,
    PERSONAS,
    PHRASES,
    absorbed,
    assistant,
    collapsed,
    expected_hash,
    head_of,
    human,
    load_path,
    make_brief_repo,
    parse_block,
    run_brief,
    tool_result,
    write_transcript,
)
from harness_fixtures import add_worktree, git

INCREMENT = "docs/increments/x.md"
STAMP = "2026-10-04T06:12:00.000Z"


@pytest.fixture
def repo(tmp_path: Path) -> Path:
    root = make_brief_repo(tmp_path.resolve() / "repo")
    write(root, INCREMENT, "# x\n\nStatus: design.\n\nNothing to quote here.\n")
    return root


@pytest.fixture
def home(tmp_path: Path, repo: Path) -> Path:
    place = tmp_path.resolve() / "home"
    write_transcript(place, repo, [human("Go ahead with the red step.", STAMP)])
    return place


@pytest.fixture
def trees(tmp_path: Path, repo: Path) -> dict[str, Path]:
    """Worktrees A, B and C of the fixture, each on its own branch.

    Each gets its own uncommitted copy of INCREMENT, as the main checkout has:
    --increment is resolved against --worktree (§3.1, §12 point 1), so a tree
    without it would be refused for that reason before the rule a test names.
    """
    base = tmp_path.resolve()
    made = {name: add_worktree(repo, base / name, f"wt-{name.lower()}") for name in "ABC"}
    for tree in made.values():
        write(tree, INCREMENT, (repo / INCREMENT).read_text())
    return made


def write(repo: Path, relative: str, text: str) -> Path:
    path = repo / relative
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    return path


def args_for(persona: str, repo: Path, *extra: str) -> list[str]:
    """The minimal valid arguments for `persona` in the main checkout."""
    base = [persona, "--worktree", str(repo), "--beside", "none"]
    if persona in NEEDS_INCREMENT:
        base += ["--increment", INCREMENT]
    return [*base, *extra]


def ok(repo: Path, home: Path, *args: str) -> str:
    result = run_brief(repo, home, *args)
    assert result.returncode == 0, f"exit {result.returncode}; stderr: {result.stderr}"
    return result.stdout


def refused(repo: Path, home: Path, *args: str, says: str) -> str:
    """Exit 2, no block, and `says` (the reason brief.py writes) in stderr.

    `says` is required: a refusal checked by exit code alone passes for any
    reason, including one raised before the rule the test names.
    """
    result = run_brief(repo, home, *args)
    assert result.returncode == 2, (
        f"exit {result.returncode}, not 2; stdout: {result.stdout}; stderr: {result.stderr}"
    )
    assert "<<<BRIEF" not in result.stdout, "a refusal printed a block"
    assert says in result.stderr, f"the refusal does not say {says!r}; stderr: {result.stderr}"
    return result.stderr


def line_starting(text: str, prefix: str) -> str:
    found = [line for line in text.splitlines() if line.startswith(prefix)]
    assert len(found) == 1, f"{len(found)} lines start with {prefix!r} in:\n{text}"
    return found[0]


# ---------------------------------------------------------------- 1. shape


@pytest.mark.parametrize("persona", PERSONAS)
def test_1_the_block_has_the_header_and_end_lines_with_one_hash(
    repo: Path, home: Path, persona: str
) -> None:
    block = parse_block(ok(repo, home, *args_for(persona, repo)))
    assert block.persona == persona
    assert block.end_hash == block.hash
    assert block.worktree == str(repo.resolve())
    assert Path(block.worktree).is_absolute()
    assert block.head == head_of(repo)


@pytest.mark.parametrize("persona", PERSONAS)
def test_1_the_hash_is_block_hash_over_the_body(repo: Path, home: Path, persona: str) -> None:
    block = parse_block(ok(repo, home, *args_for(persona, repo)))
    assert block.hash == expected_hash(persona, block.worktree, block.head, block.body)
    brief = load_path(repo / BRIEF, "h9_brief_hash")
    assert brief.block_hash(persona, block.worktree, block.head, block.body) == block.hash


def test_1_the_worktree_is_printed_absolute_and_resolved(
    repo: Path, home: Path, trees: dict[str, Path]
) -> None:
    roundabout = trees["A"] / ".." / trees["A"].name
    block = parse_block(ok(repo, home, "tester", "--worktree", str(roundabout),
                           "--beside", "none", "--increment", INCREMENT))  # fmt: skip
    assert block.worktree == str(trees["A"].resolve())
    assert block.head == head_of(trees["A"])


def test_1_two_runs_in_the_same_second_give_the_same_block(repo: Path, home: Path) -> None:
    args = args_for("tester", repo)
    for _ in range(5):
        before = datetime.now().strftime("%H%M%S")
        first, second = ok(repo, home, *args), ok(repo, home, *args)
        if datetime.now().strftime("%H%M%S") == before:
            assert first == second
            return
    pytest.fail("could not run brief.py twice within one second in five tries")


# ---------------------------------------------------------------- 2. template


@pytest.mark.parametrize("persona", PERSONAS)
def test_2_the_block_holds_the_templates_phrases_filled_in(
    repo: Path, home: Path, persona: str
) -> None:
    body = collapsed(parse_block(ok(repo, home, *args_for(persona, repo))).body)
    for phrase in PHRASES:
        assert collapsed(phrase) in body, f"{phrase!r} missing"
    assert f".claude/agents/{persona}.md from disk" in body
    assert f"You are @{persona}." in body
    assert str(repo.resolve()) in body
    assert "$" not in body, "an unfilled template variable"
    note = re.search(
        rf"{re.escape(str(repo.resolve()))}/\.claude/current-task/{persona}-\d{{6}}\.md", body
    )
    assert note is not None, "the note file is not in the block"
    if persona in NEEDS_INCREMENT:
        assert INCREMENT in body
    else:
        assert "none named" in body


def test_2_an_unknown_template_variable_is_an_error_not_a_blank(repo: Path, home: Path) -> None:
    template = repo / COMMON
    assert template.exists(), f"{COMMON} is missing from the copy"
    template.write_text(template.read_text() + "\nAlso $nonesuch here.\n")
    result = run_brief(repo, home, *args_for("tester", repo))
    assert result.returncode != 0
    assert "<<<BRIEF" not in result.stdout


# ---------------------------------------------------------------- 3. increment lines


def increment_text(body_matches: int, review: str | None) -> str:
    """Status line; `body_matches` matching body lines; an optional ## Review section."""
    lines = ["# Increment x", "", "Status: designed, round 2.", ""]
    words = ("invariant-critical", "mutation", "@perf", "acceptance run", "tools/bench.py")
    for i in range(body_matches):
        lines.append(f"Body line {i} names the {words[i % len(words)]} rule.")
    lines += ["", "An ordinary line.", "", "## Next section", "", "Plain text."]
    if review is not None:
        lines += ["", "## Review", "", review]
    return "\n".join(lines) + "\n"


REVIEW = (
    "**Round 1.** Verdict: CHANGES REQUESTED. The mutation round is missing.\n\n"
    "**Round 2.** Verdict: APPROVED. All closed."
)


def quoted(output: str) -> list[str]:
    return [line for line in output.splitlines() if re.match(r"^  \d+: ", line)]


def test_3_the_status_line_body_lines_and_review_count_are_quoted(repo: Path, home: Path) -> None:
    text = increment_text(2, REVIEW)
    write(repo, INCREMENT, text)
    output = ok(repo, home, *args_for("tester", repo))
    assert "Status: designed, round 2." in output
    numbered = {n + 1: line.strip() for n, line in enumerate(text.splitlines())}
    expected = [f"  {n}: {line}" for n, line in numbered.items() if line.startswith("Body line")]
    assert quoted(output) == expected  # the review's "mutation" line is not quoted
    assert line_starting(output, "Review rounds recorded:").startswith("Review rounds recorded: 2")
    assert "last: **Round 2.** Verdict: APPROVED. All closed." in output
    assert f"From {INCREMENT}" in output


def test_3_at_most_eight_lines_then_a_count(repo: Path, home: Path) -> None:
    write(repo, INCREMENT, increment_text(10, REVIEW))
    output = ok(repo, home, *args_for("tester", repo))
    assert len(quoted(output)) == 8
    assert "  ... 2 more; read them in the file" in output.splitlines()


def test_3_without_a_review_section_it_says_none(repo: Path, home: Path) -> None:
    write(repo, INCREMENT, increment_text(1, None))
    output = ok(repo, home, *args_for("tester", repo))
    assert "Review rounds recorded: none" in output.splitlines()


def test_3_a_long_verdict_line_is_cut_to_200_characters(repo: Path, home: Path) -> None:
    verdict = "**Round 1.** Verdict: APPROVED. " + "x" * 300
    write(repo, INCREMENT, increment_text(1, verdict))
    output = ok(repo, home, *args_for("tester", repo))
    last = output.split("last: ", 1)[1].splitlines()[0]
    assert len(last) == 200
    assert verdict.startswith(last)


# ---------------------------------------------------------------- 4. increment required


@pytest.mark.parametrize("persona", NEEDS_INCREMENT)
def test_4_an_increment_is_required_for_tester_developer_reviewer(
    repo: Path, home: Path, persona: str
) -> None:
    refused(repo, home, persona, "--worktree", str(repo), "--beside", "none",
            says=f"@{persona} needs --increment")  # fmt: skip
    refused(repo, home, persona, "--worktree", str(repo), "--beside", "none",
            "--increment", "docs/increments/missing.md",
            says="--increment docs/increments/missing.md does not exist")  # fmt: skip


def test_4_architect_may_name_an_increment_that_does_not_exist_yet(repo: Path, home: Path) -> None:
    output = ok(repo, home, "architect", "--worktree", str(repo), "--beside", "none",
                "--increment", "docs/increments/new.md")  # fmt: skip
    assert "docs/increments/new.md (new: you create it)" in output


@pytest.mark.parametrize("persona", ["perf", "orchestrator"])
def test_4_perf_and_orchestrator_need_no_increment(repo: Path, home: Path, persona: str) -> None:
    output = ok(repo, home, persona, "--worktree", str(repo), "--beside", "none")
    assert not any(line.startswith("From ") for line in output.splitlines())


def test_4_an_unknown_persona_is_refused(repo: Path, home: Path) -> None:
    refused(repo, home, "general-purpose", "--worktree", str(repo), "--beside", "none",
            says="general-purpose is not a persona")  # fmt: skip


# ---------------------------------------------------------------- 5. worktree


def test_5_a_directory_that_is_not_a_checkout_is_refused(
    repo: Path, home: Path, tmp_path: Path
) -> None:
    plain = tmp_path / "plain"
    plain.mkdir()
    refused(repo, home, "tester", "--worktree", str(plain),
            "--beside", "none", "--increment", INCREMENT,
            says=f"{plain} is not the top of a checkout")  # fmt: skip


def test_5_a_subdirectory_of_a_checkout_is_refused(repo: Path, home: Path) -> None:
    """§3.1: `rev-parse --show-toplevel` must resolve to the path itself."""
    refused(repo, home, "tester", "--worktree", str(repo / "docs"), "--beside", "none",
            "--increment", INCREMENT,
            says=f"{repo / 'docs'} is not the top of a checkout")  # fmt: skip


def test_5_a_path_with_whitespace_is_refused(repo: Path, home: Path, tmp_path: Path) -> None:
    spaced = add_worktree(repo, tmp_path.resolve() / "a b", "wt-spaced")
    refused(repo, home, "tester", "--worktree", str(spaced), "--beside", "none",
            "--increment", INCREMENT, says="contains whitespace")  # fmt: skip


def test_5_the_main_checkout_is_accepted(repo: Path, home: Path) -> None:
    assert parse_block(ok(repo, home, *args_for("tester", repo))).worktree == str(repo.resolve())


# ---------------------------------------------------------------- 6. concurrency


def concurrency(output: str) -> str:
    """The Concurrency part, its lines joined (the generated text may wrap)."""
    lines = output.splitlines()
    start = next(i for i, line in enumerate(lines) if line.startswith("Concurrency:"))
    end = next(i for i, line in enumerate(lines[start:], start) if line.startswith("Write limit:"))
    return collapsed(" ".join(lines[start:end]))


def run_as(persona: str, worktree: Path, *beside: str, extra: tuple[str, ...] = ()) -> list[str]:
    args = [persona, "--worktree", str(worktree), "--increment", INCREMENT, *extra]
    for entry in beside:
        args += ["--beside", entry]
    return args


def test_6_beside_none_runs_alone(repo: Path, home: Path) -> None:
    assert concurrency(ok(repo, home, *args_for("tester", repo))) == "Concurrency: you run alone."


def test_6_beside_is_required(repo: Path, home: Path) -> None:
    refused(repo, home, "tester", "--worktree", str(repo), "--increment", INCREMENT,
            says="the following arguments are required: --beside")  # fmt: skip


def test_6_none_with_another_entry_is_refused(
    repo: Path, home: Path, trees: dict[str, Path]
) -> None:
    refused(repo, home, *run_as("developer", trees["B"], "none", f"tester:{trees['A']}"),
            says="--beside none goes alone, not with other --beside entries")  # fmt: skip


@pytest.mark.parametrize(
    ("entry", "says"),
    [
        ("tester", "--beside tester: give <persona>:<worktree>"),
        ("tester:{missing}", "{missing} is not the top of a checkout"),
        ("nobody:{A}", "--beside nobody:{A}: give <persona>:<worktree>"),
    ],
)
def test_6_a_malformed_beside_entry_is_refused(
    repo: Path, home: Path, trees: dict[str, Path], tmp_path: Path, entry: str, says: str
) -> None:
    places = {"missing": tmp_path / "missing", "A": trees["A"]}
    text = entry.format(**places)
    refused(repo, home, *run_as("developer", trees["B"], text), says=says.format(**places))


def test_6_beside_one_writer_names_it_and_allows_a_build(
    repo: Path, home: Path, trees: dict[str, Path]
) -> None:
    line = concurrency(ok(repo, home, *run_as("developer", trees["B"], f"tester:{trees['A']}")))
    assert f"@tester (writer, in {trees['A'].resolve()})" in line
    assert "At most two writers run at once" in line
    assert "You may build C++ in this run." in line


def test_6_no_build_says_so(repo: Path, home: Path, trees: dict[str, Path]) -> None:
    args = run_as("developer", trees["B"], f"tester:{trees['A']}", extra=("--no-build",))
    line = concurrency(ok(repo, home, *args))
    assert "You may not build C++ in this run." in line
    assert "You may build" not in line


@pytest.mark.parametrize(
    ("persona", "where", "beside"),
    [
        ("reviewer", "C", ("tester:A", "developer:B")),
        ("developer", "B", ("tester:A", "reviewer:C")),
    ],
)
def test_6_a_read_only_persona_runs_beside_two_writers(
    repo: Path,
    home: Path,
    trees: dict[str, Path],
    persona: str,
    where: str,
    beside: tuple[str, ...],
) -> None:
    entries = [f"{e.split(':')[0]}:{trees[e.split(':')[1]]}" for e in beside]
    line = concurrency(ok(repo, home, *run_as(persona, trees[where], *entries)))
    if persona == "developer":
        assert f"@reviewer (read-only, in {trees['C'].resolve()})" in line
    else:
        assert f"@developer (writer, in {trees['B'].resolve()})" in line
    assert "read-only agents do not count toward the two" in line


def test_6_two_read_only_entries_beside_two_writers(
    repo: Path, home: Path, trees: dict[str, Path], tmp_path: Path
) -> None:
    d = add_worktree(repo, tmp_path.resolve() / "D", "wt-d")
    args = run_as("developer", trees["B"], f"tester:{trees['A']}",
                  f"reviewer:{trees['C']}", f"reviewer:{d}")  # fmt: skip
    line = concurrency(ok(repo, home, *args))
    assert f"@reviewer (read-only, in {d.resolve()})" in line


@pytest.mark.parametrize(
    ("persona", "beside"),
    [("perf", "tester"), ("tester", "perf"), ("perf", "reviewer"), ("reviewer", "perf")],
)
def test_6_rule_1_nothing_runs_beside_a_timing_run(
    repo: Path, home: Path, trees: dict[str, Path], persona: str, beside: str
) -> None:
    refused(repo, home, *run_as(persona, trees["B"], f"{beside}:{trees['A']}"),
            says="nothing runs beside a timing run (@perf)")  # fmt: skip


def test_6_perf_alone_says_nothing_runs_beside_it(repo: Path, home: Path) -> None:
    output = ok(repo, home, "perf", "--worktree", str(repo), "--beside", "none")
    assert concurrency(output) == ("Concurrency: you run alone; nothing runs beside a timing run.")


def test_6_rule_2_at_most_two_writers(repo: Path, home: Path, trees: dict[str, Path]) -> None:
    args = run_as("architect", trees["C"], f"tester:{trees['A']}", f"developer:{trees['B']}")
    refused(repo, home, *args, says="at most two writers run at once")


def test_6_rule_3_two_writers_in_one_worktree(
    repo: Path, home: Path, trees: dict[str, Path]
) -> None:
    refused(repo, home, *run_as("developer", trees["A"], f"tester:{trees['A']}"),
            says="two writers share one worktree")  # fmt: skip


@pytest.mark.parametrize(
    ("persona", "beside"), [("reviewer", "developer"), ("developer", "reviewer")]
)
def test_6_rule_4_a_reader_shares_no_worktree_with_a_writer(
    repo: Path, home: Path, trees: dict[str, Path], persona: str, beside: str
) -> None:
    shared = trees["A"].resolve()
    refused(repo, home, *run_as(persona, trees["A"], f"{beside}:{trees['A']}"),
            says=f"@reviewer is read-only and shares {shared} with a writer")  # fmt: skip


def test_6_read_only_is_derived_from_an_empty_write_limit(
    repo: Path, home: Path, trees: dict[str, Path], capsys: pytest.CaptureFixture[str],
    monkeypatch: pytest.MonkeyPatch,
) -> None:  # fmt: skip
    """With `orchestrator` given no write limit, it counts as read-only: beside two
    writers it is accepted and named read-only. `main(argv)` is assumed (handback)."""
    monkeypatch.setenv("HOME", str(home))
    brief = load_path(repo / BRIEF, "h9_brief_derived")
    original: Any = brief.WRITES["orchestrator"]
    monkeypatch.setitem(brief.WRITES, "orchestrator", type(original)())
    argv = ["orchestrator", "--worktree", str(trees["C"]), "--beside", f"tester:{trees['A']}",
            "--beside", f"developer:{trees['B']}"]  # fmt: skip
    try:
        code = brief.main(argv)
    except SystemExit as stop:
        code = stop.code
    out = capsys.readouterr()
    assert code in (0, None), out.err
    assert "read-only agents do not count toward the two" in collapsed(out.out)


# ---------------------------------------------------------------- 7. Ola


def with_turns(repo: Path, home: Path, *entries: dict[str, Any]) -> None:
    write_transcript(home, repo, list(entries))


def test_7_a_quotation_from_a_human_turn_is_printed_collapsed_with_its_time(
    repo: Path, home: Path
) -> None:
    with_turns(repo, home, human("Ok.  Read only agents should be allowed\n even when two agents "
                                "are working.", "2026-10-04T05:40:00.000Z"))  # fmt: skip
    output = ok(repo, home, *args_for("tester", repo),
                "--ola", "Read only agents  should be allowed even\nwhen two agents")  # fmt: skip
    assert "Ola, verbatim" in output
    assert (
        '"Read only agents should be allowed even when two agents" (2026-10-04T05:40:00.000Z)'
        in output
    )


@pytest.mark.parametrize("entry", [assistant, tool_result])
def test_7_a_quotation_found_only_outside_human_turns_is_refused(
    repo: Path, home: Path, entry: Any
) -> None:
    with_turns(
        repo, home, human("Something else.", STAMP), entry("turn off fused multiply-add", STAMP)
    )
    refused(repo, home, *args_for("tester", repo), "--ola", "turn off fused multiply-add",
            says="--ola 'turn off fused multiply-add' is not in any human turn")  # fmt: skip


def test_7_a_quotation_in_an_absorbed_queued_prompt_is_found(repo: Path, home: Path) -> None:
    with_turns(
        repo, home, human("Start.", STAMP), absorbed("yes to stage 1", "2026-10-04T06:00:00.000Z")
    )
    output = ok(repo, home, *args_for("tester", repo), "--ola", "yes to stage 1")
    assert '"yes to stage 1" (2026-10-04T06:00:00.000Z)' in output


def test_7_no_session_id_is_refused(repo: Path, home: Path) -> None:
    result = run_brief(repo, home, *args_for("tester", repo), "--ola", "Go ahead", session=None)
    assert result.returncode == 2
    assert "<<<BRIEF" not in result.stdout
    assert "--ola 'Go ahead' is not in any human turn" in result.stderr


def test_7_no_transcript_is_refused(repo: Path, tmp_path: Path) -> None:
    empty = tmp_path / "empty-home"
    empty.mkdir()
    refused(repo, empty, *args_for("tester", repo), "--ola", "Go ahead",
            says="--ola 'Go ahead' is not in any human turn")  # fmt: skip


def test_7_without_ola_there_is_no_ola_part(repo: Path, home: Path) -> None:
    output = ok(repo, home, *args_for("tester", repo))
    assert not [line for line in output.splitlines() if line.startswith("Ola, verbatim, checked")]


# ---------------------------------------------------------------- 8. note file, write limit

#: §3.1's table, one fragment per cell that must appear in WRITES[persona].
WRITE_LIMITS = {
    "tester": ("tests/",),
    "developer": ("src_python/", "include/", "src/", "bindings/", "tools/", ".claude/hooks/",
                  ".github/", "CMakeLists.txt", "pyproject.toml"),
    "perf": ("docs/benchmarks/",),
    "architect": ("docs/", "docs/retrospectives/", "ROADMAP.md", "CLAUDE.md"),
    "orchestrator": ("docs/retrospectives/",),
    "reviewer": (),
}  # fmt: skip


@pytest.mark.parametrize("persona", PERSONAS)
def test_8_the_note_file_is_named_under_the_main_checkout_and_not_created(
    repo: Path, home: Path, trees: dict[str, Path], persona: str
) -> None:
    args = [persona, "--worktree", str(trees["A"]), "--beside", "none"]
    if persona in NEEDS_INCREMENT:
        args += ["--increment", INCREMENT]
    output = ok(repo, home, *args)
    pattern = rf"{re.escape(str(repo.resolve()))}/\.claude/current-task/{persona}-(\d{{6}})\.md"
    found = re.search(pattern, output)
    assert found is not None, "the note path is not under the main checkout"
    assert not (repo / ".claude" / "current-task" / f"{persona}-{found.group(1)}.md").exists()
    assert not (trees["A"] / ".claude" / "current-task").exists()


@pytest.mark.parametrize("persona", PERSONAS)
def test_8_writes_holds_the_designs_table(repo: Path, persona: str) -> None:
    brief = load_path(repo / BRIEF, f"h9_brief_writes_{persona}")
    assert set(brief.WRITES) == set(PERSONAS)
    limit = brief.WRITES[persona]
    for fragment in WRITE_LIMITS[persona]:
        assert fragment in str(limit), f"{fragment} missing from WRITES[{persona!r}]"
    assert bool(limit) == bool(WRITE_LIMITS[persona])


@pytest.mark.parametrize("persona", PERSONAS)
def test_8_the_write_limit_line_states_writes(repo: Path, home: Path, persona: str) -> None:
    line = line_starting(ok(repo, home, *args_for(persona, repo)), "Write limit:")
    assert line.endswith("and your note file.")
    for fragment in WRITE_LIMITS[persona]:
        assert fragment in line
    if persona == "reviewer":
        others = {f for fragments in WRITE_LIMITS.values() for f in fragments}
        assert not [f for f in others if f in line], line


# ---------------------------------------------------------------- 9. hash


def test_9_any_change_to_body_persona_worktree_or_head_changes_the_hash(repo: Path) -> None:
    brief = load_path(repo / BRIEF, "h9_brief_hash9")
    head = head_of(repo)
    base = brief.block_hash("tester", str(repo), head, "line one\nline two")
    other_head = "f" * 40 if head != "f" * 40 else "e" * 40
    for changed in (
        brief.block_hash("tester", str(repo), head, "line one\nline twO"),
        brief.block_hash("developer", str(repo), head, "line one\nline two"),
        brief.block_hash("tester", str(repo / "x"), head, "line one\nline two"),
        brief.block_hash("tester", str(repo), other_head, "line one\nline two"),
    ):
        assert changed != base
    assert base == expected_hash("tester", str(repo), head, "line one\nline two")


def test_9_trailing_spaces_and_crlf_do_not_change_the_hash(repo: Path, home: Path) -> None:
    brief = load_path(repo / BRIEF, "h9_brief_hash9b")
    text = ok(repo, home, *args_for("tester", repo))
    pasted = "\r\n".join(line + "  \t" for line in text.splitlines())
    [block] = brief.find_blocks(pasted)
    assert brief.check(block, repo) is None


def test_9_one_changed_character_is_caught_by_check(repo: Path, home: Path) -> None:
    brief = load_path(repo / BRIEF, "h9_brief_hash9c")
    text = ok(repo, home, *args_for("tester", repo))
    edited = text.replace("the files win", "the files lose", 1)
    assert edited != text
    [block] = brief.find_blocks(edited)
    reason = brief.check(block, repo)
    assert reason is not None and "edited" in reason


def test_9_the_block_is_made_at_the_worktrees_head(
    repo: Path, home: Path, trees: dict[str, Path]
) -> None:
    git(trees["A"], "commit", "-q", "--allow-empty", "-m", "moved")
    block = parse_block(ok(repo, home, "tester", "--worktree", str(trees["A"]), "--beside", "none",
                           "--increment", INCREMENT))  # fmt: skip
    assert block.head == head_of(trees["A"]) != head_of(repo)


# ---------------------------------------------------- §12 point 1: the increment's checkout
#
# §3.1 as ruled after code review round 1: a relative --increment is resolved
# against --worktree, not against the checkout running brief.py; an absolute
# one is used as given; the block shows the path as given. `trees` gives each
# worktree its own copy of INCREMENT; a test that needs it absent removes it.


def test_12_the_block_quotes_the_worktrees_copy_of_the_increment(
    repo: Path, home: Path, trees: dict[str, Path]
) -> None:
    write(
        repo, INCREMENT, "# x\n\nStatus: main's status.\n\nMain's line names the mutation rule.\n"
    )
    write(
        trees["A"],
        INCREMENT,
        "# x\n\nStatus: the branch's status.\n\nThe branch's line names the @perf rule.\n",
    )
    output = ok(repo, home, "tester", "--worktree", str(trees["A"]), "--beside", "none",
                "--increment", INCREMENT)  # fmt: skip
    assert "Status: the branch's status." in output
    assert "  5: The branch's line names the @perf rule." in output.splitlines()
    assert "main's status" not in output.lower()
    assert "Main's line" not in output
    assert f"From {INCREMENT}" in output  # the path as given


def test_12_a_file_only_in_the_worktree_is_accepted_for_tester(
    repo: Path, home: Path, trees: dict[str, Path]
) -> None:
    only = "docs/increments/branch-only.md"
    write(trees["A"], only, "# y\n\nStatus: only on the branch.\n")
    assert not (repo / only).exists()
    output = ok(repo, home, "tester", "--worktree", str(trees["A"]), "--beside", "none",
                "--increment", only)  # fmt: skip
    assert "Status: only on the branch." in output


def test_12_for_architect_a_file_missing_from_the_worktree_is_new(
    repo: Path, home: Path, trees: dict[str, Path]
) -> None:
    (trees["A"] / INCREMENT).unlink()
    assert (repo / INCREMENT).exists()
    output = ok(repo, home, "architect", "--worktree", str(trees["A"]), "--beside", "none",
                "--increment", INCREMENT)  # fmt: skip
    assert f"{INCREMENT} (new: you create it)" in output


def test_12_a_file_in_the_main_checkout_only_is_refused_for_tester(
    repo: Path, home: Path, trees: dict[str, Path]
) -> None:
    """The converse of the branch-only case: master's copy does not stand in."""
    (trees["A"] / INCREMENT).unlink()
    assert (repo / INCREMENT).exists()
    refused(repo, home, "tester", "--worktree", str(trees["A"]), "--beside", "none",
            "--increment", INCREMENT, says=f"--increment {INCREMENT} does not exist")  # fmt: skip


def test_12_an_absolute_increment_is_used_as_given(
    repo: Path, home: Path, trees: dict[str, Path], tmp_path: Path
) -> None:
    elsewhere = write(
        tmp_path.resolve() / "elsewhere", "inc.md", "# z\n\nStatus: from elsewhere.\n"
    )
    write(trees["A"], "inc.md", "# z\n\nStatus: the worktree's inc.\n")
    output = ok(repo, home, "tester", "--worktree", str(trees["A"]), "--beside", "none",
                "--increment", str(elsewhere))  # fmt: skip
    assert "Status: from elsewhere." in output
    assert "the worktree's inc" not in output
    assert f"From {elsewhere}" in output
