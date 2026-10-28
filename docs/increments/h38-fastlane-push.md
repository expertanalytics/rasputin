# h38: fast-lane push and enqueue without a fresh yes, present or away

Status: design, awaiting Ola's approval (§8). No production code yet.

Ola, verbatim, 2026-10-28, on whether "present" covers him asleep without
`away.py`: "I set away when going to sleep."

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
- Two holes in away mode today, found by @reviewer in round 1 and probed
  again in a scratch repository (a bare `origin`, a clone, branch `quick`
  adding `src/x.cpp`, a full-lane path), calling today's code in process:
  - `fast_lane()` diffs against the *local* ref `refs/remotes/origin/master`.
    Before: `False`. After writing the branch's head into
    `.git/refs/remotes/origin/master` with `printf` (which no guard stops):
    `True`, while `git ls-remote origin refs/heads/master` still named the
    real master. So today a full-lane branch passes the away-mode fast lane.
  - `home()` decides "holds private commits" from the local
    `refs/remotes/private/master`. A branch built on a private-only commit:
    `('private', '')`. After writing an unrelated one-commit history into
    that ref: `('origin', '')`. So today a granted branch holding private
    commits may be pushed to origin. `review_of()`, `recorded()` and
    `changes()` read the same local refs (a forged one there lets a merge of
    an unreviewed branch pass as a merge of master).
- `.github/` is neither a full-lane path nor governed, so a fast-lane branch
  may change what CI runs, `CI result` included. #367, #369 and #374 changed
  `main.yaml` last week.

## 3. The design: one path, five changes to it

All in `tools/push_grant.py`, on the existing `fast_lane()` / `judge()` path,
plus the off switch in `tools/away.py`. No second path, no new module.

```
guard_push.py --(command, reasons)--> judge
  judge: mode off | on (broken: None, as now)
    flag  = unattended.json only when on; {} when off      (change 1)
    real  = masters(where): ls-remote of each remote       (change 4)
    act   = granted_act(command)                            (unchanged)
    fast  = not granted and fast_lane(where, state, act, real)
              no full-lane, increment, governed or CI-config path (change 2)
              against real["origin"], <= 700 lines
    granted path: home / review_of / recorded / changes, against real
    create: real head of the branch on origin == local head (change 4)
    enqueue: checks(pr_view(where, branch, repo), head)      (change 3)
  -> Grant passed (logged in grants-used.jsonl) | refused (asks, or queued away) | None (asks)
```

1. **Present too.** `judge` no longer returns `None` when away mode is off:
   with mode `off`, only the fast-lane branch is judged, and the flag is
   taken as `{}` (no grants, no window), so `unattended.json` is read only
   when the mode is `on` (round 1: it was read unconditionally and would
   raise when absent, which the hook turns into a block). A present
   enqueue is logged with `"window": null`, and "once per head" then means
   once per head while present. With mode `broken` it still returns `None`.
   "Present" is what the guard can see: away mode off. Settled: Ola sets
   away mode when he goes to sleep (his words, top of this file).
2. **No governed or CI-config path.** `fast_lane()` also fails when any path
   in the diff is one `guard_governance.governed()` names (imported from
   `.claude/hooks/`: rule files, `.claude/settings*.json*`, agents, hooks,
   skills, `tools/check_*`, and the tools CI runs that are not `check_*`:
   `ci_changes.py`, `count_loc.py`; `main.yaml` runs no other `tools/` script,
   checked with `grep -on "tools/[a-z_]*\.py" .github/workflows/main.yaml`),
   or one in a new tuple `CI_CONFIG` in `push_grant.py`: the prefix `.github/`,
   and any file whose name is one CI's tools read as their settings,
   anywhere in the tree: `pyproject.toml`, `conftest.py`, `CMakeLists.txt`,
   `CMakePresets.json`, and the files that would take precedence over
   `pyproject.toml` if added: `pytest.ini`, `.pytest.ini`, `tox.ini`,
   `setup.cfg`, `mypy.ini`, `.mypy.ini`, `ruff.toml`, `.ruff.toml`,
   `.coveragerc`. Such a branch asks as before. Not in the list: test files.
   A test module can still skip itself, but only a `conftest.py` or a
   setting can make every suite pass; a fast-lane writer writes its own
   tests (`docs/increments/README.md`, *Two lanes*).
3. **CI green before enqueue.** A fast-lane `enqueue` passes only when
   `checks(pr_view(where, branch, repo))` finds the pull request open into
   master, its head equal to the branch head, and every check green with at
   least one success. `pr_view` takes the repository as an argument
   (round 1: it was hard-wired to the private one): `PUBLIC_REPO =
   "expertanalytics/rasputin"` for origin, `PRIVATE_REPO` as now for the
   private merge. Otherwise the refusal reason is appended and the act asks
   (present) or is refused and queued (away).
4. **The remotes' real master, not the local refs.** A new
   `masters(where) -> (dict[remote, sha], reason)` runs, once per judged
   act, `git ls-remote <url> refs/heads/master` for origin, and for
   `private` when `remote.private.url` is set, where `<url>` is
   `git remote get-url --push <remote>` (which expands `insteadOf` and
   `pushInsteadOf`, so the guard reads the repository the push would write).
   For a create it also asks origin for `refs/heads/<branch>` in the same
   call. A master sha missing locally is fetched once by the guard
   (`git fetch -q <remote> refs/heads/master`); still missing, a timeout
   (15 s per call) or any git error refuses, with the reason. Every place
   `push_grant.py` reads `refs/remotes/<remote>/master` today (`fast_lane`,
   the `count_loc.py` base, `home`, `changes`, `review_of`, `recorded`)
   takes that sha instead. This fixes both holes of §2 on both paths, the
   granted one included.
   Also the round-1 follow-up: `gh pr create --head <b>` opens the PR from
   origin's `<b>`, not the local branch the guard judged. A create now
   passes only when origin's `refs/heads/<b>` equals the local head.
   Cost, measured: `git ls-remote` of origin and of the private repository
   each took 1.7 to 1.9 s (three runs, from this machine); so a judged act
   costs about 3.6 s more, on a guard that already runs `count_loc.py`. Only
   push and `gh pr` commands reach `judge`.
5. **Off switch** (waiting on Ola, question 3; included by default).
   `python3 tools/away.py --fast-push off` writes
   `<git-common-dir>/harness/fast-push.json`; while that file exists (or
   cannot be read: fails closed) `fast_lane()` returns false in both modes,
   so every fast-lane act asks again. `away.py --fast-push on` deletes it,
   after a `y` typed on the terminal, as entering away mode needs (it
   loosens). `off` needs no terminal, so `!python3 tools/away.py --fast-push
   off` works at a session prompt. Agents cannot run `away.py` or write
   `harness/` (`guard_governance`, unchanged). Default with no file: on.

Changes 2 to 4 apply in away mode too: one path, one set of conditions
(question 1).

**Acts let through**, on origin only, for a branch passing all of the above:
the push of its head by sha; `gh pr create` of it into master (waiting on
Ola, question 2; included by default); and one `gh pr merge` enqueue per head.
The writer waits for CI with `gh pr checks <pr> --watch` run in the
background, then enqueues; an early enqueue just asks.

**What stays a question:** every full-lane or increment-file branch; any
governed or CI-config path; over 700 lines; the private remote; a force
push or history rewrite; any command not in the three bare forms (an `&&`
chain, a variable); `gh pr edit|ready|update-branch`; an enqueue with CI
red, pending or unreadable, or a second enqueue of the same head; a stale
named sha; a create whose origin branch is not the local head; a remote
that cannot be read; any `.claude/settings*.json` edit (`guard_governance`,
unchanged); everything while the switch is off.

**Left open, on purpose.** A global `~/.gitconfig` can point `insteadOf`
for reading at one repository and `pushInsteadOf` for writing at another;
`get-url --push` then names the second, which is the one written, so the
judgement is about the right repository. Not closed: a test module that
disables others from import time (no file list can find it).

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

If Ola answers no to question 2 (create still asks), the added sentence is
instead "A fast-lane branch for origin is pushed and enqueued with no yes,
grant or review: the guard checks lane and CI." (21 words, net 0). Question
3 changes no rule text: the off switch and the path lists live in the
guard's code, since no agent acts on them; the guard does.

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
`guard_governance.governed()` names, no `CI_CONFIG` path, and
`count_loc.py --max 700` passing. Script (scratch, not committed), re-run
in round 2 with the `CI_CONFIG` class added.

| class (PRs) | wait, sum | in waits over 1 h |
|---|---|---|
| fast lane (35: retrospectives and harness fixes) | 101.4 h | 94.3 h in 15 PRs |
| CI config, would now ask (3: #367, #369, #374) | 8.6 h | 7.6 h in 1 |
| governed path, would still ask (9) | 9.0 h | 5.6 h |
| full lane or increment file, untouched (28) | 54.9 h | 48.6 h |

Most long fast-lane waits span a night (#386 gil-flake 11.2 h); the guard
passes those in either mode.
What remains after the change is CI on the PR before the enqueue: 0.2 h
median (17 PRs opened at once, at most 0.6 h). Subtracting 0.3 h per PR
leaves **about 90 h over the week, 2.6 h per fast-lane PR, an upper bound**
(a wait also holds time Ola spent reading, which the change does not
remove). Since a retrospective counts as an item in flight (`CLAUDE.md` §3,
three at most), most of that time was also a writer slot held open.

What it does not touch: the mesh and land-surface PRs are full lane (#339
ar50-vegetation 17.3 h, #326 network-catchment 3.5 h); their 54.9 h still
waits on Ola's yes.

## 6. Tests (`tests/python/test_push_grant.py`, `@tester`)

Fixture change first. `masters()` asks the remotes themselves, so every
fixture that now sets `refs/remotes/<remote>/master` with `update-ref` gets
a real bare repository on a local path as that remote, and pushes its
master there (which also sets the tracking ref): `harness_fixtures.py`
(the sites at lines 429, 455, 499, 502, 530 today), `test_push_grant.py`
(284, 704, and `Fast` at 1200), `test_away.py` (1073), `test_guard_push.py`
(1195). Measured: `git ls-remote` of a local bare repository took 0.10 s.
`Fast` becomes: `master` (1 commit, from `make_repo`) pushed to a bare
`origin`, branch `quick` with one commit of `files`, pushed to origin's
`quick` too unless a test says otherwise. With `tools/y.py` (`y = 1`), the
diff is 1 path and 1 net line. Each test judges in process; "present" is
the flag file removed.

1. Present, fast branch: push and create pass; each logs one
   `grants-used.jsonl` line with `"lane": "fast"` and `"window": null`.
   Replaces `test_a_fast_lane_branch_goes_nowhere_while_ola_is_present`.
2. Present, enqueue with a stubbed `pr_view` reporting green checks on the
   head: passes, and the stub saw `expertanalytics/rasputin`; a second
   enqueue of the same head is refused.
3. Present and away, enqueue with checks pending, failed, on another head,
   none reported, or unreadable: refused, reason names it.
4. Present and away, a branch adding each of `CLAUDE.md`,
   `.claude/settings.json`, `.claude/agents/x.md`, `.claude/hooks/x.py`,
   `tools/check_x.py`, `tools/push_grant.py`, `.github/workflows/main.yaml`,
   `pyproject.toml`, `tests/python/conftest.py`, `tests/cpp/CMakeLists.txt`,
   `pytest.ini`, `ruff.toml`: not judged (`None`), nothing logged. Amends
   `test_a_fast_lane_branch_passes_without_grant_or_review` to drop
   `CLAUDE.md` from its files.
5. Present, the existing not-fast cases (full lane, increment file, 701
   lines, another sha): still `None`.
6. Forged refs, both modes: (a) branch adding `src/x.cpp`, then
   `refs/remotes/origin/master` written to its head through the file:
   `None`. (b) Away, granted branch built on a private-only commit, then
   `refs/remotes/private/master` written to an unrelated one-commit
   history: a push to origin is refused, reason "holds commits of
   private/master". (c) Origin's master moved in the bare repository and
   not fetched: the guard fetches it and judges against it (a branch whose
   diff against the new master is fast passes).
7. Create when origin's `quick` is absent or at another sha: refused, the
   reason names origin's sha or its absence.
8. Origin's URL a missing path: every act `None` or refused, nothing logged.
9. Switch file present (and, separately, a directory in its place): every
   fast act `None`, both modes.
10. Mode broken: `None`.
11. `away.py --fast-push off` writes the file with no terminal; `on` without
    a `y` on the terminal leaves it; with `y` removes it.
12. `guard_push` end to end, present, no `unattended.json` on disk: the fast
    push line exits 0 with no output; a full-lane push still asks.

Planting check: each of 1, 3, 4, 6a, 6b, 7 and 12 is run once against
`push_grant.py` with its new condition removed (for 12: the unconditional
read put back), and must fail.

## 7. Size

Production: `push_grant.py` about 50 lines (`masters` about 20, `CI_CONFIG`
and its check 8, the threaded sha and the create check about 12, present
mode and `pr_view` 10), `away.py` about 20, `guard_push.py` 1: about 70 net
lines, all in governed files, so the build needs Ola's yes for each write
(as any harness change). Tests about 180 lines, fixture changes included.
Lane: full (it changes the approval boundary), `@reviewer` before push.

## 8. Questions for Ola

1. Today, while away, a fast-lane branch that changes a rule file such as
   `CLAUDE.md` may be pushed and enqueued, and its enqueue does not wait
   for CI. Apply the new conditions to away mode too? Default: yes.
2. Opening the pull request (`gh pr create`) was not in your words "push
   and enqueue". Let it through without a yes too? Default: yes (in
   `tools/questions.py`, item fastlane-push; waiting on the answer).
3. An off switch you can flip from a session prompt (`!python3
   tools/away.py --fast-push off`), which you did not ask for. Build it?
   Default: yes (in `tools/questions.py`, item fastlane-push; waiting on
   the answer).

## Review

- Design round 1: CHANGES REQUESTED, 22d09fe1..140f9e11, 0 (bypasses tried and closed: full-lane path, increment file, over 700 lines, stale sha, force push, other branch or master, private remote or `-R`, governed file, CI not green, enqueue for another head; blocking: `.github/` is neither full lane nor governed, so a fast-lane branch can make `CI result` always green (#367, #369, #374 changed `main.yaml`): exclude `.github/` and maybe `pyproject.toml` in `fast_lane()`, add to test 4; a forged `refs/remotes/origin/master` (Write or `printf >` passes both guards) empties the lane diff, and a forged `private/master` can send private commits to origin: check against `git ls-remote origin refs/heads/master`; the `away.py --fast-push` switch and `gh pr create` go beyond Ola's words; `judge` reads `unattended.json` unconditionally and would crash when he is present; `pr_view` is hard-wired to the private repository; follow-up: `gh pr create` is judged on the local branch, the PR opens from the remote one; measurement re-run: 38 PRs, 109.9 h, an upper bound; rule text -2 words holds).
