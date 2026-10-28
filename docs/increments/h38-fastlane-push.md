# h38: fast-lane push and enqueue without a fresh yes, present or away

Status: design, awaiting Ola's approval (§8). No production code yet.

Ola, 2026-10-28, on the proposal "Fast-lane push and enqueue without asking
while you are present, as away mode already allows: CI green, lane check green,
and never for rule files or settings": "Decisions 1-4: defaults."

## 1. Prior art: legacy and literature

*Literature.* This is the "not rocket science rule" (Graydon Hoare, 2014,
<https://graydon2.dreamwidth.org/1597.html>): a change reaches the main branch
only through an automated test of the result, never by a person's say-so.
GitHub's merge queue (already ours, h10) implements it. The one departure:
we also put the decision to *publish* (push, open a pull request) behind a
mechanical check, where those tools leave it to the author. That is
stricter than the method, not looser. No novelty is claimed.

*Legacy.* `git grep -n -i -l "merge queue\|auto-merge\|git push"
legacy-archive -- legacy` returns nothing; nothing is carried across.

## 2. What exists (read, not assumed)

- `.claude/hooks/guard_push.py` asks before every push and every
  `gh pr create|merge`; before asking it calls `push_grant.judge`, and a
  passed judgement lets the act through silently.
- `tools/push_grant.py` `judge` returns `None` (so: ask) unless away mode is
  on. While on, `fast_lane()` passes a branch for origin whose diff against
  `origin/master` touches no `check_lane.FULL_LANE` path and no increment
  file, at most 700 net lines (`count_loc.py --max 700`), named in one of the
  three bare forms `granted_act` accepts:
  `git push origin <sha>:refs/heads/<b>` (the head's sha, no other option, so
  never a force push), `gh pr create --head <b> ...` (base master), and
  `gh pr merge <b> --match-head-commit <sha>` (once per head per window,
  logged in `<git-common-dir>/harness/grants-used.jsonl`).
- Two gaps against Ola's wording, both in away mode today: a fast-lane branch
  that changes a rule file passes (the suite pins it:
  `test_a_fast_lane_branch_passes_without_grant_or_review` changes
  `CLAUDE.md`), and a fast-lane enqueue does not look at CI (the merge queue
  is relied on).

## 3. The design: one path, three changes to it

All in `tools/push_grant.py`, on the existing `fast_lane()` / `judge()` path.
No second path, no new module.

1. **Present too.** `judge` no longer returns `None` when away mode is off:
   with mode `off`, only the fast-lane branch is judged (no grants exist);
   with mode `broken` it still returns `None` (refused, as now). "Present"
   is what the guard can see: away mode off. That includes Ola asleep without
   `away.py`; together with the away rule this makes the fast lane pass
   at any time.
2. **No governed path.** `fast_lane()` also fails when any path in the diff
   is one `guard_governance.governed()` names (imported from
   `.claude/hooks/`, the one list: rule files, `.claude/settings*.json*`,
   agents, hooks, skills, `tools/check_*`, the self-protecting tools such as
   `push_grant.py` itself). Such a branch asks as before.
3. **CI green before enqueue.** A fast-lane `enqueue` passes only when
   `checks(pr_view(...))` (the private-merge check, given the repository as
   an argument) finds the pull request open into master, its head equal to
   the branch head, and every check green with at least one success.
   Otherwise the refusal reason is appended and the act asks (present) or is
   refused and queued (away), as any refusal is now.

Changes 2 and 3 apply in away mode too: one path, one set of conditions
(question 1 for Ola).

**Acts let through**, on origin only, for a branch passing all of the above:
the push of its head by sha, `gh pr create` of it into master, and one
`gh pr merge` enqueue per head.
The writer waits for CI with `gh pr checks <pr> --watch` run in the
background, then enqueues; an early enqueue just asks.

**What stays a question:** every full-lane or increment-file branch; any
governed path; over 700 lines; the private remote; a force push or history
rewrite; any command not in the three bare forms (an `&&` chain, a variable);
`gh pr edit|ready|update-branch`; an enqueue with CI red, pending or
unreadable, or a second enqueue of the same head; a stale named sha; any
`.claude/settings*.json` edit (`guard_governance`, unchanged); everything
while the switch below is off.

**How Ola turns it off.** `python3 tools/away.py --fast-push off` writes
`<git-common-dir>/harness/fast-push.json`; while that file exists (or cannot
be read: fails closed) `fast_lane()` returns false in both modes, so every
fast-lane act asks again. `away.py --fast-push on` deletes it, after a `y`
typed on the terminal, as entering away mode needs (it loosens). `off` needs
no terminal, so `!python3 tools/away.py --fast-push off` works at a session
prompt. Agents cannot run `away.py` or write `harness/` (`guard_governance`,
unchanged). Default with no file: on.

`guard_push.py` changes one word: its refusal says "fast-lane branch" rather
than "granted branch" when the lane was fast.

## 4. Rule text (`.claude/REQUIRED-READING.md`), net words -2

- In **A fresh yes, every time**, after "...or the permission system.", add:
  "A fast-lane branch for origin (`docs/increments/README.md`) needs no yes,
  grant or review; the guard checks its lane and CI." (+19 words, `wc -w`)
- In **Unattended mode**, delete "A fast-lane branch for origin
  (`docs/increments/README.md`) needs neither grant nor review; the guard
  checks its lane." (-16)
- In `docs/increments/README.md`, **Two lanes**, delete "The approval
  boundary (`.claude/REQUIRED-READING.md`) holds." (-5): REQUIRED-READING is
  the one statement of it.

The off switch and the governed-path rule live in the guard's code, not in
rule text: no agent acts on them, the guard does.

## 5. What it would have saved

Source: `python3 docs/retrospectives/lead_time.py 2026-10-21` (75 PRs merged
10-21 to 10-28). Its "wait" runs from the branch's last commit to the merge
commit's committer time, which on a merge-queue repository is the moment the
PR entered the queue (#394: queue entry 19:30 UTC, merge commit 19:30 UTC),
so the queue's own CI run is outside it and stays either way.

How fast-lane PRs were identified: the fast lane as a rule began only on
10-28 (#389), so each PR was classified after the fact by the guard's own
test as this design leaves it: `git diff --name-only <merge>^1...<merge>^2`
with no `check_lane.FULL_LANE` path, no increment file, no path
`guard_governance.governed()` names, and `count_loc.py --max 700` passing.
Script (scratch, not committed): it also read each PR's opening time and its
first `added_to_merge_queue` event from GitHub.

| class (PRs) | wait, sum | in waits over 1 h |
|---|---|---|
| fast lane (38: 32 retrospectives, 6 harness/CI fixes) | 110.0 h | 101.9 h in 16 PRs |
| governed path, would still ask (9) | 9.0 h | 5.6 h |
| full lane or increment file, untouched (28) | 54.9 h | 48.6 h |

Of the 16 long fast-lane waits, 11 span a night (#386 gil-flake 11.2 h,
#374 ci-queue-serial 7.6 h); the guard passes those too, since it sees only
whether away mode is on (question 2). What remains after the change is CI on the PR before
the enqueue: 0.2 h median (17 PRs opened at once, at most 0.6 h). Subtracting
0.3 h per PR leaves **about 100 h saved over the week, 2.7 h per fast-lane
PR**. Since a retrospective counts as an item in flight (`CLAUDE.md` §3,
three at most), most of that time was also a writer slot held open.

What it does not touch: the mesh and land-surface PRs are full lane (#339
ar50-vegetation 17.3 h, #326 network-catchment 3.5 h); their 54.9 h still waits on Ola's yes.

## 6. Tests (`tests/python/test_push_grant.py`, `@tester`)

Each judges in process on the existing `Fast` fixture (a two-commit repository,
`origin/master` plus branch `quick`), the flag file removed for "present".

1. Present, fast branch (`tools/y.py`, 1 line): push and create pass; each logs
   one `grants-used.jsonl` line with `"lane": "fast"` and no window.
   Replaces `test_a_fast_lane_branch_goes_nowhere_while_ola_is_present`.
2. Present, enqueue with a stubbed `pr_view` reporting green checks on the
   head: passes; a second enqueue of the same head is refused (asks).
3. Present and away, enqueue with checks pending, failed, on another head,
   none reported, or unreadable: refused, reason names it.
4. Present and away, a branch adding each of `CLAUDE.md`,
   `.claude/settings.json`, `.claude/agents/x.md`, `.claude/hooks/x.py`,
   `tools/check_x.py`, `tools/push_grant.py`: not judged (`None`), nothing
   logged. Amends `test_a_fast_lane_branch_passes_without_grant_or_review`
   to drop `CLAUDE.md` from its files.
5. Present, the existing not-fast cases (full lane, increment file, 701 lines,
   another sha): still `None`.
6. Switch file present (and, separately, a directory in its place): every
   fast act `None`, both modes.
7. Mode broken: `None`.
8. `away.py --fast-push off` writes the file with no terminal; `on` without
   a `y` on the terminal leaves it; with `y` removes it.
9. `guard_push` end to end, present: the fast push line exits 0 with no
   output; a full-lane push still asks.

Planting check: each of 1, 3, 4, 6 is run once against `push_grant.py`
with its new condition removed, and must fail.

## 7. Size

Production: `push_grant.py` about 25 lines, `away.py` about 20,
`guard_push.py` 1: about 45 net lines, all in governed files, so the build
needs Ola's yes for each write (as any harness change). Tests about 130
lines. Lane: full (it changes the approval boundary), `@reviewer` before push.

## 8. Questions for Ola

1. Today, while away, a fast-lane branch that changes a rule file such as
   `CLAUDE.md` may be pushed and enqueued, and its enqueue does not wait
   for CI. Apply the two new conditions to away mode too? Default: yes.
2. "Present" can only mean "away mode off", so a fast-lane branch also
   goes out while you sleep without `away.py`. Accept? Default: yes.

## Review

- Design round 1: CHANGES REQUESTED, 22d09fe1..140f9e11, 0 (bypasses tried and closed: full-lane path, increment file, over 700 lines, stale sha, force push, other branch or master, private remote or `-R`, governed file, CI not green, enqueue for another head; blocking: `.github/` is neither full lane nor governed, so a fast-lane branch can make `CI result` always green (#367, #369, #374 changed `main.yaml`): exclude `.github/` and maybe `pyproject.toml` in `fast_lane()`, add to test 4; a forged `refs/remotes/origin/master` (Write or `printf >` passes both guards) empties the lane diff, and a forged `private/master` can send private commits to origin: check against `git ls-remote origin refs/heads/master`; the `away.py --fast-push` switch and `gh pr create` go beyond Ola's words; `judge` reads `unattended.json` unconditionally and would crash when he is present; `pr_view` is hard-wired to the private repository; follow-up: `gh pr create` is judged on the local branch, the PR opens from the remote one; measurement re-run: 38 PRs, 109.9 h, an upper bound; rule text -2 words holds).
