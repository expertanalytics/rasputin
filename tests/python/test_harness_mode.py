"""`tools/harness_mode.py`: the mode, `decide`, and the queue (h3 U1, T1-T3).

The spec is `docs/increments/h3-unattended-u1.md`: §3.2 (the flag and
`read_mode`), §3.4 (`decide` and the queue entry), §4 (T1-T3). The mode only
ever restricts: an unreadable flag is `broken`, and `broken` refuses.
"""

from __future__ import annotations

import json
import os
from datetime import UTC, datetime, timedelta
from pathlib import Path
from typing import Any

import pytest

from harness_fixtures import (
    GUARD_PUSH,
    MAX_WINDOW,
    Tool,
    bash_event,
    flag_broken,
    flag_on,
    flag_path,
    flag_window,
    iso,
    make_repo,
    pretool_decision,
    queue_lines,
    run_script,
    state_dir,
    write_flag,
)

# h17 §4b: a harness test; CI runs it in the `harness` job, not the product legs.
pytestmark = pytest.mark.harness

NOW = datetime(2026, 9, 30, 18, 0, 0, tzinfo=UTC)
UNTIL = datetime(2026, 10, 1, 3, 36, 0, tzinfo=UTC)
ACT = "git push origin feature"
WHY = "git push writes to the remote"
GUARD_REASON = "This reaches beyond the working tree: git push writes to the remote."

needs_permissions = pytest.mark.skipif(
    hasattr(os, "geteuid") and os.geteuid() == 0, reason="root ignores file permissions"
)


hm = Tool("harness_mode")


@pytest.fixture
def repo(tmp_path: Path) -> Path:
    return make_repo(tmp_path / "repo")


@pytest.fixture
def chmodded() -> Any:
    touched: list[Path] = []
    yield touched
    for path in touched:
        if path.exists():
            path.chmod(0o700)


# ---------------------------------------------------------------- T1 read_mode


def test_state_dir_is_the_git_common_dir_plus_harness(repo: Path) -> None:
    found = hm.state_dir(repo)
    assert found is not None
    assert Path(found).resolve() == state_dir(repo).resolve()


def test_state_dir_outside_a_repository_is_none(tmp_path: Path) -> None:
    lonely = tmp_path / "not-a-repo"
    lonely.mkdir()
    assert hm.state_dir(lonely) is None


def test_read_mode_without_a_state_dir_is_broken() -> None:
    mode = hm.read_mode(None, NOW)
    assert mode.state == "broken"
    assert mode.detail == "git common dir not found"


def test_read_mode_with_no_flag_is_off(repo: Path) -> None:
    assert hm.read_mode(state_dir(repo), NOW).state == "off"


def test_read_mode_with_a_valid_window_is_on_with_its_until(repo: Path) -> None:
    flag_window(repo, NOW, UNTIL)
    mode = hm.read_mode(state_dir(repo), NOW + timedelta(hours=1))
    assert mode.state == "on"
    assert mode.until == UNTIL


def test_read_mode_ignores_set_by_and_keep_awake_pid(repo: Path) -> None:
    write_flag(repo, {"since": iso(NOW), "until": iso(UNTIL), "set_by": 7, "keep_awake_pid": "x"})
    assert hm.read_mode(state_dir(repo), NOW).state == "on"


def test_read_mode_at_exactly_until_is_off(repo: Path) -> None:
    flag_window(repo, NOW, UNTIL)
    assert hm.read_mode(state_dir(repo), UNTIL).state == "off"


def test_read_mode_after_until_is_off_and_leaves_the_file(repo: Path) -> None:
    flag_window(repo, NOW, UNTIL)
    assert hm.read_mode(state_dir(repo), UNTIL + timedelta(seconds=1)).state == "off"
    assert flag_path(repo).exists()


def test_read_mode_accepts_a_window_of_exactly_the_cap_plus_60_s(repo: Path) -> None:
    flag_window(repo, NOW, NOW + MAX_WINDOW + timedelta(seconds=60))
    assert hm.read_mode(state_dir(repo), NOW).state == "on"


BROKEN_CONTENT: dict[str, str | dict[str, Any]] = {
    "not JSON": "{not json",
    "a JSON list": '["since", "until"]',
    "since missing": {"until": iso(UNTIL)},
    "until missing": {"since": iso(NOW)},
    "since not a string": {"since": 1727719200, "until": iso(UNTIL)},
    "until not parseable": {"since": iso(NOW), "until": "tomorrow morning"},
    "naive times": {
        "since": iso(NOW.replace(tzinfo=None)),
        "until": iso(UNTIL.replace(tzinfo=None)),
    },
    "until equal to since": {"since": iso(NOW), "until": iso(NOW)},
    "until before since": {"since": iso(UNTIL), "until": iso(NOW)},
    "window over the cap": {
        "since": iso(NOW),
        "until": iso(NOW + MAX_WINDOW + timedelta(seconds=61)),
    },
}


@pytest.mark.parametrize("content", BROKEN_CONTENT.values(), ids=BROKEN_CONTENT.keys())
def test_read_mode_on_an_invalid_flag_is_broken_with_a_detail(
    repo: Path, content: str | dict[str, Any]
) -> None:
    write_flag(repo, content)
    mode = hm.read_mode(state_dir(repo), NOW)
    assert mode.state == "broken"
    assert isinstance(mode.detail, str) and mode.detail.strip()


@pytest.mark.parametrize("kind", ["directory", "naive", "reversed", "too-long", "not-json"])
def test_read_mode_on_each_broken_fixture_flag_is_broken(repo: Path, kind: str) -> None:
    flag_broken(repo, kind, now=NOW)
    mode = hm.read_mode(state_dir(repo), NOW)
    assert mode.state == "broken"
    assert mode.detail.strip()


@needs_permissions
def test_read_mode_on_an_unreadable_flag_is_broken(repo: Path, chmodded: list[Path]) -> None:
    path = flag_window(repo, NOW, UNTIL)
    path.chmod(0o000)
    chmodded.append(path)
    mode = hm.read_mode(state_dir(repo), NOW)
    assert mode.state == "broken"
    assert mode.detail.strip()


# ---------------------------------------------------------------- T2 decide


def _mode(state: str) -> Any:
    if state == "on":
        return hm.Mode(state="on", until=UNTIL, detail="")
    if state == "broken":
        return hm.Mode(state="broken", until=None, detail="not a JSON object")
    return hm.Mode(state="off", until=None, detail="")


@pytest.mark.parametrize("state", ["off", "on", "broken"])
def test_decide_pass_passes_in_every_mode(state: str) -> None:
    decision = hm.decide("pass", "", _mode(state), ACT, WHY)
    assert decision.kind == "pass"
    assert decision.queue is False


@pytest.mark.parametrize("state", ["off", "on", "broken"])
def test_decide_deny_keeps_the_guards_reason_and_never_queues(state: str) -> None:
    decision = hm.decide("deny", "not an agent's act", _mode(state), ACT, WHY)
    assert (decision.kind, decision.reason, decision.queue) == (
        "deny",
        "not an agent's act",
        False,
    )


def test_decide_ask_while_attended_asks_with_the_guards_reason() -> None:
    decision = hm.decide("ask", GUARD_REASON, _mode("off"), ACT, WHY)
    assert (decision.kind, decision.reason, decision.queue) == ("ask", GUARD_REASON, False)


def test_decide_ask_while_unattended_denies_with_the_queue_reason() -> None:
    decision = hm.decide("ask", GUARD_REASON, _mode("on"), ACT, WHY)
    assert decision.kind == "deny"
    assert decision.queue is True
    assert decision.reason.startswith(
        f"Refused: unattended mode is on until 2026-10-01 03:36 UTC; {ACT} waits for Ola ({WHY})."
    )
    for phrase in (
        "The refusal is already recorded in the queue.",
        "Do not retry",
        "ASK OLA:",
        "handback",
    ):
        assert phrase in decision.reason


def test_decide_ask_with_a_broken_flag_denies_with_the_broken_reason() -> None:
    decision = hm.decide("ask", GUARD_REASON, _mode("broken"), ACT, WHY)
    assert decision.kind == "deny"
    assert decision.queue is True
    assert decision.reason.startswith(
        "Refused: the unattended flag cannot be read (not a JSON object), so every guarded "
        "act is refused until Ola runs python3 tools/away.py --back; "
        f"{ACT} waits for Ola ({WHY})."
    )
    assert "ASK OLA:" in decision.reason
    assert "unattended mode is on until" not in decision.reason


# ---------------------------------------------------------------- T3 the queue


ENTRY_FIELDS = {"at", "branch", "cwd", "agent_type", "agent_id", "hook", "act", "why"}


def test_queue_entry_in_the_main_session_has_null_agent_fields(repo: Path) -> None:
    entry = hm.queue_entry(bash_event(repo, ACT), "guard_push", ACT, WHY, NOW)
    assert set(entry) >= ENTRY_FIELDS
    assert datetime.fromisoformat(entry["at"]) == NOW
    assert entry["branch"] == "master"
    assert entry["cwd"] == str(repo)
    assert (entry["agent_type"], entry["agent_id"]) == (None, None)
    assert (entry["hook"], entry["act"], entry["why"]) == ("guard_push", ACT, WHY)


def test_queue_entry_in_a_subagent_carries_its_type_and_id(repo: Path) -> None:
    event = bash_event(repo, ACT, agent_id="a1", agent_type="tester")
    entry = hm.queue_entry(event, "guard_push", ACT, WHY, NOW)
    assert (entry["agent_type"], entry["agent_id"]) == ("tester", "a1")


def test_queue_entry_outside_a_repository_has_a_null_branch(tmp_path: Path) -> None:
    entry = hm.queue_entry(bash_event(tmp_path, ACT), "guard_push", ACT, WHY, NOW)
    assert entry["branch"] is None


def test_queue_entry_truncates_the_act_to_1000_characters(repo: Path) -> None:
    long_act = "echo " + "x" * 3000
    entry = hm.queue_entry(bash_event(repo, long_act), "guard_push", long_act, WHY, NOW)
    assert len(entry["act"]) <= 1000
    assert entry["act"] == long_act[: len(entry["act"])]


def test_append_queue_creates_the_dir_and_writes_one_line_per_call(repo: Path) -> None:
    first = hm.queue_entry(bash_event(repo, ACT), "guard_push", ACT, WHY, NOW)
    second = hm.queue_entry(
        bash_event(repo, "gh pr merge 1"),
        "guard_push",
        "gh pr merge 1",
        "gh pr changes a pull request",
        NOW,
    )
    assert not state_dir(repo).exists()
    assert hm.append_queue(state_dir(repo), first) is None
    assert hm.append_queue(state_dir(repo), second) is None
    assert queue_lines(repo) == [json.loads(json.dumps(first)), json.loads(json.dumps(second))]
    assert (state_dir(repo) / "queue.jsonl").stat().st_mode & 0o777 == 0o600


def test_append_queue_without_a_state_dir_returns_an_error() -> None:
    entry = {"act": ACT}
    error = hm.append_queue(None, entry)
    assert isinstance(error, str) and error.strip()


@needs_permissions
def test_append_queue_into_a_read_only_dir_returns_an_error(
    repo: Path, chmodded: list[Path]
) -> None:
    state_dir(repo).mkdir()
    state_dir(repo).chmod(0o500)
    chmodded.append(state_dir(repo))
    entry = hm.queue_entry(bash_event(repo, ACT), "guard_push", ACT, WHY, NOW)
    error = hm.append_queue(state_dir(repo), entry)
    assert isinstance(error, str) and error.strip()


def test_two_refusals_through_the_hook_give_two_queue_lines(repo: Path) -> None:
    flag_on(repo)
    for command in ("git push", "gh pr merge 1"):
        found = pretool_decision(run_script(repo, GUARD_PUSH, bash_event(repo, command)))
        assert found is not None and found[0] == "deny"
    lines = queue_lines(repo)
    assert [line["act"] for line in lines] == ["git push", "gh pr merge 1"]
    for line in lines:
        assert set(line) >= ENTRY_FIELDS
        assert line["hook"] == "guard_push"
        assert line["branch"] == "master"


@needs_permissions
def test_a_refusal_the_queue_cannot_record_is_still_a_deny_and_says_so(
    repo: Path, chmodded: list[Path]
) -> None:
    flag_on(repo)
    state_dir(repo).chmod(0o500)
    chmodded.append(state_dir(repo))
    found = pretool_decision(run_script(repo, GUARD_PUSH, bash_event(repo, "git push")))
    assert found is not None
    kind, reason = found
    assert kind == "deny"
    assert "could not be recorded" in reason
    assert "your ASK OLA line is the only record" in reason
    assert "The refusal is already recorded in the queue." not in reason
