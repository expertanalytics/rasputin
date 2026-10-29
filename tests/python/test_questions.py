"""`tools/questions.py`: the questions for Ola, checked, numbered and answered.

The spec is `docs/increments/h35-question-tool.md`, *4. Tests, `tests/python/test_questions.py`*.
Every test runs `command()` on a store in `tmp_path`, so the live store under
the git common directory is never read or written; the clock is the fixed
`NOW` passed in. Inputs are literal strings.
"""

from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import pytest

import harness_fixtures

# A harness test; CI runs it in the `harness` job, not the product legs.
pytestmark = pytest.mark.harness

questions = harness_fixtures.Tool("questions")

NOW = datetime(2026, 10, 27, 18, 0, tzinfo=UTC)


# ---------------------------------------------------------------- helpers


@pytest.fixture
def store(tmp_path: Path) -> Path:
    """The store path; its directory holds everything the tool writes."""
    return tmp_path / "common" / "questions.jsonl"


def run(store: Path, *argv: str) -> list[str]:
    """The lines one command prints."""
    return list(questions.command(questions.parser().parse_args(list(argv)), store, NOW))


def refused(store: Path, *argv: str) -> str:
    """The refusal's sentence; the store is left as it was."""
    before = snapshot(store.parent)
    with pytest.raises(questions.RefusedError) as caught:
        run(store, *argv)
    assert snapshot(store.parent) == before
    return str(caught.value)


def snapshot(directory: Path) -> dict[str, bytes]:
    if not directory.is_dir():
        return {}
    return {str(p.relative_to(directory)): p.read_bytes()
            for p in sorted(directory.rglob("*")) if p.is_file()}  # fmt: skip


def records(store: Path) -> list[dict[str, Any]]:
    good, bad = questions.load(store)
    assert bad == 0
    return list(good)


def newest(store: Path, text: str) -> dict[str, Any]:
    """The newest record whose text is `text`."""
    return [r for r in records(store) if r["text"] == text][-1]


def numbered(lines: list[str]) -> list[str]:
    return [line for line in lines if line[:1].isdigit()]


def ask(store: Path, text: str, *extra: str, default: str = "no") -> list[str]:
    return run(store, "add", text, "--default", default, *extra)


# ---------------------------------------------------------------- 1. C1 default


@pytest.mark.parametrize("default", ["", "TBD", "none", "?", "-"])
def test_c1_a_question_without_a_real_default_is_refused(store: Path, default: str) -> None:
    reason = refused(store, "add", "Cache the repair?", "--default", default)
    assert "needs a default" in reason


def test_c1_a_real_default_is_accepted(store: Path) -> None:
    ask(store, "Cache the repair?", default="no")
    assert newest(store, "Cache the repair?")["default"] == "no"


def test_c1_taken_needs_a_default_too(store: Path) -> None:
    assert "needs a default" in refused(store, "taken", "the cache stays local", "--default", "")


# ---------------------------------------------------------------- 2. C2 one decision


def test_c2_two_question_marks_are_refused(store: Path) -> None:
    assert "one decision" in refused(store, "add", "Push it? And merge?", "--default", "no")


@pytest.mark.parametrize("extra", [[], ["--act"]])
def test_c2_a_question_must_end_with_a_question_mark(store: Path, extra: list[str]) -> None:
    reason = refused(store, "add", "Push the tool.", "--default", "no", *extra)
    assert "end with" in reason


def test_c2_one_question_mark_at_the_end_is_accepted(store: Path) -> None:
    ask(store, "Push the tool?")
    assert newest(store, "Push the tool?")["state"] == "open"


# ---------------------------------------------------------------- 3. no question is held


def test_a_second_yes_no_question_is_open_and_listed_with_the_first(store: Path) -> None:
    ask(store, "Cache the repair?")
    said = ask(store, "Push the tool and open its pull request?", "--act")
    assert said == ["stored as q2"]
    ask(store, "Which unit for heights?", "--choices", "metre|foot", default="metre")
    ask(store, "Which unit for depths?", "--choices", "metre|fathom", default="metre")
    assert {r["state"] for r in records(store)} == {"open"}
    shown = run(store, "list")
    assert shown == numbered(shown) and len(shown) == 4, shown
    assert not any("held" in line for line in shown)
    assert run(store, "answer", "4", "yes") == ['Ola, 2026-10-27: "yes"']
    assert run(store, "drop", "1", "superseded") == ["Dropped: Cache the repair?"]
    assert newest(store, "Push the tool and open its pull request?")["state"] == "answered"


def test_a_stored_held_record_is_listed_and_answered_as_open(store: Path) -> None:
    ask(store, "Cache the repair?")
    ask(store, "Push the tool?")
    held = {**newest(store, "Push the tool?"), "state": "held"}
    with store.open("a") as out:
        out.write(json.dumps(held) + "\n")
    shown = run(store, "list")
    assert [line[:2] for line in shown] == ["1.", "2."], shown
    assert shown[1].startswith("2. Push the tool?")
    refused(store, "add", "Push the tool?", "--default", "no")
    run(store, "answer", "2", "yes")
    assert newest(store, "Push the tool?")["state"] == "answered"


# ---------------------------------------------------------------- 4. C4 self-contained


@pytest.mark.parametrize("text", ["Same as above?", "As in question 2?"])
def test_c4_a_question_pointing_elsewhere_is_refused(store: Path, text: str) -> None:
    assert "read alone" in refused(store, "add", text, "--default", "no")


@pytest.mark.parametrize("text", ["Keep cells above 1000 m?", "Is it the same as master?"])
def test_c4_bare_above_and_same_as_pass(store: Path, text: str) -> None:
    ask(store, text)
    assert newest(store, text)["state"] == "open"


def test_c4_phrases_match_only_as_whole_words(store: Path) -> None:
    ask(store, "Keep what Ola has before ruled?")
    assert newest(store, "Keep what Ola has before ruled?")["state"] == "open"
    assert "read alone" in refused(store, "add", "Keep the order as before?", "--default", "no")


# ---------------------------------------------------------------- 5. C5 refusals


@pytest.mark.parametrize(
    ("text", "token"),
    [
        ("Ship the h33 change now?", "h33"),
        ("Do the P4 step first?", "P4"),
        ("Merge #381 now?", "#381"),
        ("Revert 7c2aa7c0 now?", "7c2aa7c0"),
        ("Ask @someone to look?", "@someone"),
        ("Push worktree-x today?", "worktree-x"),
        ("Split tools/questions.py in two?", "tools/questions.py"),
        ("Move the notes under docs/ now?", "docs/"),
        ("Keep a/b.py as it is?", "a/b.py"),
        ("Use the ABCD format?", "ABCD"),
    ],
)  # fmt: skip
def test_c5_a_label_or_internal_name_is_refused_by_name(store: Path, text: str, token: str) -> None:
    assert token in refused(store, "add", text, "--default", "no")


def test_c5_a_branch_name_is_allowed_on_an_act(store: Path) -> None:
    ask(store, "Push worktree-x today?", "--act", default="wait for your yes")
    record = newest(store, "Push worktree-x today?")
    assert (record["kind"], record["state"]) == ("act", "open")


# ---------------------------------------------------------------- 6. C5 acceptances


@pytest.mark.parametrize(
    "text",
    [
        "Build with -O1 for the timing run?",
        "Cap the mesh at 1400000 triangles?",
        "Let @perf time the change?",
        "Write the list to /Users/x/tools/questions.py?",
        "Read the settings from ~/x/y.toml?",
        "Read the sources from $RASPUTIN_DATA/sources.toml?",
        "Fetch https://a.b/c.py before the run?",
        "Report wind in m/s?",
        "Is a yes/no answer enough here?",
        "Show radiation in W/m²?",
        "Keep the switch on/off by default?",
        "Add one line to CLAUDE.md about it?",
        "Use Large Shapes (LS) for the demo, with LS as the default?",
    ],
)  # fmt: skip
def test_c5_plain_tokens_that_look_like_labels_pass(store: Path, text: str) -> None:
    ask(store, text)
    assert newest(store, text)["state"] == "open"


@pytest.mark.parametrize("text", ["Let @perf's run decide the timing?",
                                  "Let @perf\u2019s run decide the timing?"])  # fmt: skip
def test_c5_a_possessive_persona_handle_passes(store: Path, text: str) -> None:
    ask(store, text)
    assert newest(store, text)["state"] == "open"


def test_c5_a_proper_name_passes_only_when_named_with_term(store: Path) -> None:
    text = "Fetch the rail lines from Bane NOR?"
    assert "NOR" in refused(store, "add", text, "--default", "no")
    ask(store, text, "--term", "NOR")
    assert newest(store, text)["state"] == "open"


# ---------------------------------------------------------------- 7. C6 short


def sized(total: int, default: str = "no") -> str:
    """A plain question whose length plus the default's is `total` characters."""
    head, tail = "Keep the ", "?"
    return head + "x" * (total - len(default) - len(head) - len(tail)) + tail


def test_c6_240_characters_with_the_default_are_accepted(store: Path) -> None:
    text = sized(240)
    assert len(text) + len("no") == 240
    ask(store, text)
    assert newest(store, text)["state"] == "open"


def test_c6_241_characters_are_refused(store: Path) -> None:
    assert "at most 240 characters" in refused(store, "add", sized(241), "--default", "no")


# ---------------------------------------------------------------- 8. C7 not re-asked


def test_c7_re_asking_an_answered_question_quotes_the_answer(store: Path) -> None:
    ask(store, "Cache the repair?")
    run(store, "list")
    run(store, "answer", "1", "yes, cache it per area")
    reason = refused(store, "add", "  cache   the REPAIR?", "--default", "no")
    assert "yes, cache it per area" in reason
    assert "2026-10-27" in reason


def test_c7_re_asking_an_open_question_is_refused(store: Path) -> None:
    ask(store, "Cache the repair?")
    refused(store, "add", "Cache the repair?", "--default", "yes")


# ---------------------------------------------------------------- 9. list


def test_list_numbers_questions_acts_then_taken_in_one_sequence(store: Path) -> None:
    ask(store, "Push the tool and open its pull request?", "--act", default="wait for your yes")
    ask(store, "Which unit for heights?", "--choices", "metre|foot", default="metre")
    run(store, "taken", "the cache lives in the user's cache folder", "--default", "yes")
    ask(store, "Cache the repair?")
    ask(store, "Repair on every run?")
    shown = run(store, "list")
    assert shown == numbered(shown) and len(shown) == 5, shown
    assert shown[0].startswith("1. Which unit for heights?") and "Default: metre" in shown[0]
    assert shown[1].startswith("2. Cache the repair?")
    assert shown[2].startswith("3. Repair on every run?")
    assert shown[3].startswith("4. Push the tool and open its pull request?")
    assert "Default: wait for your yes" in shown[3]
    assert shown[4].startswith("5. Taken: ")
    assert "the cache lives in the user's cache folder" in shown[4]


def test_a_taken_line_appears_in_one_list_only(store: Path) -> None:
    run(store, "taken", "the cache lives in the user's cache folder", "--default", "yes")
    assert any("Taken: " in line for line in run(store, "list"))
    assert not any("Taken: " in line for line in run(store, "list"))


# ---------------------------------------------------------------- 10. list_lines is read-only


def test_list_lines_between_list_and_answer_changes_nothing(store: Path) -> None:
    ask(store, "Push the tool and open its pull request?", "--act", default="wait for your yes")
    first = run(store, "list")
    assert first[0].startswith("1. Push the tool")
    ask(store, "Which unit for heights?", "--choices", "metre|foot", default="metre")
    run(store, "taken", "the cache lives in the user's cache folder", "--default", "yes")
    before = snapshot(store.parent)
    recap = questions.list_lines(questions.load(store)[0])
    assert snapshot(store.parent) == before
    assert recap[0].startswith("1. Which unit for heights?")
    assert any("Taken: " in line for line in recap)
    run(store, "answer", "1", "yes")
    assert newest(store, "Push the tool and open its pull request?")["state"] == "answered"
    assert newest(store, "Which unit for heights?")["state"] == "open"
    after = run(store, "list")
    assert any("Taken: the cache lives" in line for line in after)


# ---------------------------------------------------------------- 11. answer


def test_answer_resolves_through_the_last_list_after_a_later_add(store: Path) -> None:
    ask(store, "Which unit for heights?", "--choices", "metre|foot", default="metre")
    run(store, "list")
    ask(store, "Which colour for rivers?", "--choices", "blue|green", default="blue")
    run(store, "answer", "1", "foot")
    assert newest(store, "Which unit for heights?")["state"] == "answered"
    assert newest(store, "Which colour for rivers?")["state"] == "open"


def test_answer_refuses_a_number_the_last_list_did_not_print(store: Path) -> None:
    ask(store, "Cache the repair?")
    refused(store, "answer", "1", "yes")
    run(store, "list")
    refused(store, "answer", "2", "yes")


def test_answer_stores_the_words_verbatim(store: Path) -> None:
    words = 'yes, but "only" for Norway\nand Sweden later'
    ask(store, "Cache the repair?")
    run(store, "list")
    printed = run(store, "answer", "1", words)
    assert f'Ola, 2026-10-27: "{words}"' in printed
    record = newest(store, "Cache the repair?")
    assert (record["state"], record["words"]) == ("answered", words)


# ---------------------------------------------------------------- 12. drop


def test_drop_closes_the_numbered_question(store: Path) -> None:
    ask(store, "Cache the repair?")
    run(store, "list")
    run(store, "drop", "1", "superseded by the per-area cache")
    assert newest(store, "Cache the repair?")["state"] == "dropped"
    assert not any("Cache the repair?" in line for line in run(store, "list"))


# ---------------------------------------------------------------- 13. the store


def test_ids_count_up_and_the_newest_record_per_id_wins(store: Path) -> None:
    ask(store, "Which unit for heights?", "--choices", "metre|foot", default="metre")
    ask(store, "Which colour for rivers?", "--choices", "blue|green", default="blue")
    assert [r["id"] for r in records(store)] == ["q1", "q2"]
    run(store, "list")
    run(store, "answer", "1", "foot")
    assert numbered(run(store, "list"))[0].startswith("1. Which colour for rivers?")


def test_a_malformed_line_is_skipped_and_counted(store: Path) -> None:
    ask(store, "Cache the repair?")
    with store.open("a") as out:
        out.write("{not json\n")
    good, bad = questions.load(store)
    assert bad == 1
    assert [r["text"] for r in good] == ["Cache the repair?"]
    assert json.loads(store.read_text().splitlines()[0])["id"] == "q1"
    assert numbered(run(store, "list"))[0].startswith("1. Cache the repair?")
