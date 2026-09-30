# Unattended mode: design

Status: design only (@architect, 2026-09-30). No code, hook, setting or rule
file is changed by this document. It extends the generic layer designed in
`docs/research/generic-harness.md` (the profile, §2-§3; R-B and the ledger,
§6) and slots into that document's §7 migration plan (§8 below). Claims
about Claude Code are cited to its documentation pages, fetched on
2026-09-30: `hooks`, `sub-agents` and `permission-modes` under
`https://code.claude.com/docs/en/`.

**Decision in one paragraph.** Unattended mode is a flag with an expiry
time. Only Ola can set it, with a command that needs a real terminal. While
it is on, every guard that would ask **denies with a reason** instead. The
hook itself records the refused act in a queue, and the reason tells the
agent to write an `ASK OLA:` line and carry on. On a non-default branch
outside the live checkout, a rule-file edit is **allowed and logged**
instead of asked. What keeps that safe is a push gate that **does not
depend on the mode**. Whenever the ref being pushed or merged changes a
rule file relative to the base branch, `guard_push` refuses until Ola has
acknowledged that exact diff. The diff is identified by the file's git blob
at the merge-base and at the head, and Ola acknowledges it with a review
command that also needs a real terminal. Leaving the mode never loosens
anything, so an agent gains nothing by leaving it. Entering the mode
loosens one thing, and that loosening is closed again at push time.

---

## 1. What the platform does with each decision (checked)

| Fact | Source |
|---|---|
| A PreToolUse `"deny"` prevents the call, and its `permissionDecisionReason` is **shown to Claude**. For `"ask"` the reason is shown to the user but **not** to Claude | `hooks`, *PreToolUse decision control* |
| Several hooks, several decisions: `deny` > `defer` > `ask` > `allow` | same |
| A hook's `"ask"` forces a prompt even in auto mode (since v2.1.211); the classifier can deny but not approve silently | same |
| `"defer"` works only with `-p`; "in interactive sessions it logs a warning and ignores the hook result" | `hooks`, *Defer a tool call for later* |
| Foreground subagent: "Permission prompts are passed through to you". Background subagent: the prompt surfaces in the main session and that subagent waits | `sub-agents`, *Run subagents in foreground or background* |
| Inside a subagent the hook input carries `agent_id` (only there) and `agent_type`; every event carries `cwd` and `permission_mode` | `hooks`, *Common input fields* |
| `PermissionRequest` runs when Claude Code is about to show a permission prompt, and can `deny` with a `message` to Claude | `hooks`, *PermissionRequest* |
| `SubagentStart` cannot block, but it can inject `additionalContext` into the subagent | `hooks`, *SubagentStart* |
| `AskUserQuestion` is never auto-approved in any mode | `permission-modes`, *Actions no mode auto-approves* |
| Writes under `.claude/` (except `.claude/worktrees`) and `.git` are protected paths: they prompt in `default` mode, go to the classifier in `auto` mode, and are denied in `dontAsk` mode | `permission-modes`, *Protected paths* |
| The Bash tool has no terminal: `sys.stdin.isatty()` is `False` and `open('/dev/tty')` fails with `ENXIO` | run on this branch from the Bash tool; in Ola's own terminal both succeed |

What follows from these facts:

- A deny inside a subagent reaches the subagent as a tool error with the
  reason attached. The subagent keeps running and can hand the refusal back.
  An ask inside a subagent stalls, whether the subagent runs in the
  foreground or the background. So in unattended mode no guard may ever
  return `ask`.
- `defer` cannot be the mechanism, because these sessions are interactive.
- A `dontAsk` launch would stop the stalls, but it denies every Bash command
  that is not allowlisted, it can only be chosen at launch, and it turns a
  hook's `ask` into a denial whose reason Claude never sees. That is why the
  guards have to deny themselves, with reasons of their own.

---

## 2. Entering and leaving

**The flag.** `<git-common-dir>/harness/unattended.json`, for example
`/Users/skavhaug/projects/rasputin/.git/harness/unattended.json`. It is
untracked by construction, and every worktree of the repository shares it
(`git rev-parse --git-common-dir` gives the same path from any of them). It
holds `{"since": <UTC ISO>, "until": <UTC ISO>, "set_by": "away", "keep_awake_pid": <int|null>}`.
The mode is on if and only if the file parses and `now < until`.

| Direction | Who | How | Why this direction is safe |
|---|---|---|---|
| **Enter** | **Ola only** | `python3 tools/away.py <duration>` in his own terminal. The command reads a `y` from `/dev/tty`, writes the flag with `until = now + duration × buffer` (`profile.unattended.buffer`, 1.2 in rasputin, which matches the keep-awake routine), and starts the keep-awake command for the same length. `profile.unattended.max_hours` caps the duration, and the duration is mandatory | Entering loosens one thing: rule-file edits on a branch are no longer asked (§3). An agent entering the mode to skip that prompt would be making a self-authored rule change. So entering must be Ola's |
| **Leave early** | Ola | `python3 tools/away.py --back`. It deletes the flag, stops the keep-awake process it started, and prints the morning summary (§5) | — |
| **Leave on expiry** | nobody | The hooks read `until` on every call | — |
| **Agent writes the flag** | never | Every guard **denies** any tool write that names `harness/` under the git dir, and any Bash command that runs `away.py` or the review command. This holds in both modes. The `/dev/tty` read is the second barrier | Leaving would dodge nothing: the push gate (§4) is the same in both modes, and in attended mode an `ask` waits for a human. A leave would only turn tonight's refusals into stalls. The deny exists so that the flag has **one writer**, not because leaving is dangerous |

**Keep-awake.** The keep-awake command (`profile.unattended.keep_awake`,
`caffeinate -is -t {seconds}` on macOS) and the mode are separate
mechanisms that one command starts together. The main session may still
start keep-awake on its own, as it does today when Ola leaves without
saying for how long. It may never set the mode. When Ola says he is
leaving, the main session asks how long, then asks him to run
`tools/away.py <duration>` in a terminal. (That routine lives in Ola's user
memory, which the generic layer does not own; §9, Q2.)

**The tripwire caveat.** An agent could get a terminal with `script`,
`expect` or a pty module. The guards scan for those next to `away.py` and
the review command, but the scan reads the command as text, so like the
other guards it is a tripwire and not a sandbox.

---

## 3. Each hook in the mode

One shared function, `harness_mode.decide(event, verdict, reason)`, maps a
guard's verdict to the decision it returns. Each guard computes a verdict
of `pass`, `ask` or `deny`, and never has to know which mode it is in:

- **attended, main session**: `ask` stays `ask`;
- **attended, subagent** (`agent_id` present): `ask` stays `ask` today, and
  becomes `deny` once R-B lands (as generic-harness §6.1 decides);
- **unattended, anywhere**: `ask` becomes `deny`. The hook appends the
  refused act to the queue (§5) and returns the queue reason below.

**The queue reason** is the same for every guard; `{act}` and `{why}` are
filled in by the guard:

> Refused: unattended mode is on until {until}; {act} waits for Ola ({why}).
> The refusal is already recorded in the queue. Do not retry it and do not
> work around it. Main session: add one line `ASK OLA: {act}: <what you
> wanted and why>` to `.claude/current-task/session.md`, then continue with
> other work. Subagent: put that `ASK OLA:` line in your handback and carry
> on with the rest of your brief. If the brief cannot proceed without this
> act, hand back now.

| Hook | Attended (today, unchanged) | Unattended |
|---|---|---|
| `guard_push`: push, `gh pr create/merge/ready/edit`, `gh release`, `gh repo …` | **ack gate first** (§4): `deny` if any rule-file diff on the ref is unacknowledged or rejected; otherwise `ask` | ack gate, then `deny` + queue. Pushes and merges always wait for Ola |
| `guard_push`: `--no-verify`, `rebase`, `reset --hard`, `filter-branch`, `commit --amend` | `ask` | `deny` + queue. The agent makes a new commit instead (§9, Q7) |
| `guard_governance`, **self-protecting set**: `.claude/settings*.json*`, `.claude/hooks/`, the harness state under the git dir, `tools/away.py`, the review command, `.claude/profile.toml` | `ask` (state dir, `away.py` and review command: `deny`) | **`deny`** + queue, on every branch. These files take effect in the running session before any review could see them |
| `guard_governance`, **other rule files**, target on the default branch, a detached HEAD, or inside the checkout the hooks run from (`$CLAUDE_PROJECT_DIR`) | `ask` | `deny` + queue |
| `guard_governance`, **other rule files**, target on any other branch in another worktree | `ask`. After the edit, a PostToolUse arm records a *prompt ack* (§4) | **pass**, returning nothing, so the normal permission flow runs, plus `additionalContext`: "Logged as pending. It binds nobody until Ola acknowledges it: do not cite or follow it." The hook appends to the pending journal (§4) |
| `guard_governance`, Bash arm (substring scan) | as today: it cannot parse shell, so a path built by variable expansion slips past it, and a read with a redirect trips it | as today, with two additions. (a) A write-shaped command that names an absolute path outside the cwd's git toplevel is `deny`, because the branch check cannot see which tree it targets. (b) A miss costs only the journal entry: the push gate derives the diff from git, not from the journal |
| `gates_after_commit` (PostToolUse) | reports; exits 2 when red | unchanged |
| `SessionStart` (recap) | recap | recap, headed by `UNATTENDED until <until>`, the queue rules above, the queue, and the pending and rejected rule-file diffs per branch |
| **new** `PreToolUse` on `AskUserQuestion` | pass | `deny` + queue ("write the question as an `ASK OLA:` line") |
| **new** `PermissionRequest` catch-all | no hook output (normal prompt) | `deny` with the queue reason. This catches the prompts no guard issues: permission rules, a protected-path write the classifier sent back, a first read outside the working directories |
| **new** `SubagentStart` | nothing | `additionalContext`: the mode, the queue rule, and the unacknowledged rule-file changes on this branch, marked "not in force" |

**Composition with R-B.** The per-persona guard (generic-harness §6.1)
produces verdicts in the same way and goes through the same `decide`.
Claude Code resolves several hooks by precedence, `deny` > `ask`, so R-B's
subagent `deny` wins over an unattended `pass`. For example, `@developer`
may not edit `CLAUDE.md` in either mode. R-B's main-session `ask` on
production code becomes `deny` + queue while unattended.

**Worktrees.** The settings run the hooks from `$CLAUDE_PROJECT_DIR` (the
main checkout), but every tree-dependent question is answered from the
event. For Edit and Write, the branch and toplevel come from `file_path`'s
own git toplevel. For Bash they come from `cwd` (rule (a) above covers the
rest). The flag, the queue, the journal and the acks all live in the shared
git common dir, so they are the same from every worktree.

**Failure direction.** If the flag exists but cannot be parsed, or the
common dir cannot be found, every guarded act is denied, and the reason
tells Ola to run `away.py --back`. Deny never stalls and never publishes.

---
## 4. The pending journal, the acknowledgement record, and the push gate

Everything lives in `<git-common-dir>/harness/`. It is untracked (§9, Q6),
shared by every worktree, and written only by hooks and by Ola's two
commands. Tool writes to it are denied (§2).

| File | Written by | One JSON line per | Fields |
|---|---|---|---|
| `pending.jsonl` (the journal) | `guard_governance`, when an unattended rule-file edit passes | attempted edit | `at`, `branch`, `worktree`, `path`, `session_id`, `agent_type`, `tool`, `command` or `file_path` |
| `acks.jsonl` (the record) | the review command (`by: "review"`); the attended PostToolUse arm (`by: "prompt"`) | decision on one diff | `path`, `base_blob`, `head_blob`, `decision` (`ack` / `reject`), `reason`, `by`, `at`, `branch`, `head_commit` |
| `queue.jsonl` | every guard that turns `ask` into `deny` | refused act | `at`, `branch`, `cwd`, `agent_type`, `hook`, `act` |

**Identity of a diff.** A rule-file change is the pair
`(base_blob, head_blob)` for one path. `base_blob` is
`git rev-parse <merge-base>:<path>` and `head_blob` is
`git rev-parse <head>:<path>`. A file that is absent on one side is written
as the null id. Git already addresses file content by hash, so no separate
hash is needed. The branch name is recorded for information only: an ack
covers the diff wherever it appears, and renaming a branch loses nothing.

**The authoritative set is derived from git, never read from the journal.**
For a ref *R* and base *B* (`profile.base_branch`, compared against the
remote-tracking ref), the set is:

```
changed(R) = { p in git diff --name-only $(git merge-base B R) R
               : governed(p) under the manifest at the merge-base OR at R }
```

The manifest is taken from both sides, because a branch that removes a
file from the governed set has made a rule change of its own. `governed`
is the one function `guard_governance` also uses. After step 4 of
generic-harness §7 it reads the manifest from `profile.rule_files` plus the
generic list, and `.claude/profile.toml` must itself be in that list.

**The gate**, `harness_acks.check(R)`, is run by `guard_push` before every
push, `gh pr create`, `gh pr ready` and `gh pr merge`, in both modes. For
each path in `changed(R)` it applies these rules, and the first one that
matches refuses:

| Condition | Refusal |
|---|---|
| the file is modified in the working tree or the index | "commit or discard it first; only committed content can be acknowledged" |
| its blob at the merge-base differs from its blob at `B` now | "the base moved under this rule file: merge `B` in, then review again" |
| the newest record for `(path, base_blob, head_blob)` is a `reject` | the reason Ola gave, and "rework or revert this file" |
| there is no `ack` record for `(path, base_blob, head_blob)` | "run the review command" |

For `gh pr merge` the hook resolves *R* to the PR's head
(`gh pr view --json headRefOid`) and refuses if that is not the local
branch head, because Ola reviewed the local branch. **A file that changes
after it was acknowledged** gets a new `head_blob`, so the old ack does not
match it and the gate refuses. **A revert** to the base content removes the
path from `changed(R)`, so nothing is left to acknowledge.

**The attended prompt ack.** In attended mode, when `guard_governance`
returned `ask` for an Edit or Write, it stores the `tool_use_id`. When
PostToolUse fires for that id, the edit ran, which means Ola approved it.
The hook then records `by: "prompt"` with the file's current
working-tree blob (`git hash-object`) and the current merge-base blob. If
the file is committed as it stood, the ack matches. If a later edit slips
past the tripwire, the blob differs and the gate asks for a review. A
PostToolUse with no stored id never records an ack, and neither does any
call made while the mode was on at PreToolUse time.

---

## 5. The review command and the morning

`python3 tools/harness_review.py [branch …]`, in Ola's terminal. With no
arguments it covers every local branch whose `changed()` set is not empty.
`away.py --back` prints the summary below and then offers to run it.

1. **Summary.** The queue, grouped by branch (what was refused, from which
   persona). The `ASK OLA:` lines from `session.md`. The rule-file diffs per
   branch with their status: unacknowledged, acknowledged, rejected, or
   stale (the base moved).
2. **Per file, one at a time:** the header (path, branch, `base_blob..head_blob`,
   how many commits touch it); the journal entries for that file (when,
   which persona, which tool); **"no journal entry" in bold** when git shows
   a change the tripwire never logged; then `git diff <merge-base> <head> -- <path>`
   through the pager.
3. **The decision is read from `/dev/tty`:** `a` acknowledges, `r` rejects
   and asks for a one-line reason, `s` skips, `q` quits. Each `a` or `r` is
   appended to `acks.jsonl` at once, so quitting midway loses nothing. The
   command refuses a file that has uncommitted changes.
4. **After the files:** a line per branch saying whether
   `harness_acks.check` now passes. It publishes nothing. A push is still a
   separate yes, given at `guard_push`'s prompt.

**A rejection sends the branch back.** The recap prints each `reject` as
`REWORK <branch> <path>: <reason>`, next to the `ASK OLA` lines. The next
session treats it as an instruction from Ola. The rework is a new commit,
so a new `head_blob` and a new review, or a revert, which leaves nothing to
review. `queue.jsonl` is archived to `queue-<date>.jsonl` by `--back` once
it has been printed.

---
## 6. Failure modes

| Failure | What happens | Why it is contained |
|---|---|---|
| **The flag is left on** while Ola is back | His next push is refused with "unattended until …; run `away.py --back`". Rule-file edits on branches stop asking | Deny never publishes. The unasked edits still meet the push gate. The recap's first line says `UNATTENDED`. The expiry puts an upper bound on the mistake |
| **The flag expires** while Ola is still away | The guards go back to `ask`, and the next guarded act stalls until he returns | This is the failure mode the mode was built to avoid, but it fails safe. The buffer and `max_hours` size it. The recap and the `SubagentStart` context state `until`, so the session can plan its work around it |
| **The clock** | Only the expiry reads the clock, in UTC. Sleep, timezones and travel to an exhibition do not affect it. Moving the clock needs root | Nothing that grants anything depends on time. Acks are bound to blobs, not timestamps |
| **The Mac sleeps** | The session freezes; nothing is lost | keep-awake is started by the same command. `-s` holds on AC only, and `away.py` says so when `pmset -g batt` reports battery |
| **A merge in the morning with a stale ack** | Two cases: `B` moved and changed the same rule file, or the branch changed the file after the ack. The gate's second or fourth row refuses | The merge-base check means Ola never merges a rule diff computed against a base he did not see |
| **A rule-file edit hidden in an unrelated branch**, or made through a path the Bash tripwire could not read | `changed(R)` comes from `git diff`, so the gate catches it. The review shows it with "no journal entry", which also says that the tripwire missed it | The gate does not depend on the journal, on the mode, or on the branch's stated purpose |
| **An edit to the live checkout** (`$CLAUDE_PROJECT_DIR`) or to the self-protecting set | Denied and queued in unattended mode | These change tonight's enforcement or tonight's subagents before any review |
| **An agent treats a pending rule as in force** | It is told otherwise twice: the edit's `additionalContext` and `SubagentStart` both mark it "not in force" | A persona that reads rules through its pointer (generic-harness §5c) would otherwise follow a rule nobody has authored |
| **The queue write races** with a parallel hook | One JSON line per `O_APPEND` write, under `PIPE_BUF` | The rule of one agent at a time makes a race unlikely anyway |
| **The agent obtains a pty** and runs `away.py` or the review command | The Bash tripwire denies `script`, `expect`, `unbuffer` and `pty` near those names, and does so in both modes | This is a tripwire, like every other guard here. The boundary is still the agent's to keep |

---

## 7. The generic/profile split

| Generic (harness repo) | Profile (`profile.toml` / copier answers) | rasputin value |
|---|---|---|
| The mode's semantics; `harness_mode.decide`; the queue reason | `owner`, rendered into the reason and the `ASK <OWNER>:` marker | Ola |
| The flag's location and format; `away.py`; the `/dev/tty` confirmation | `unattended.max_hours`, `unattended.buffer`, `unattended.keep_awake` | 72, 1.2, `caffeinate -is -t {seconds}` |
| The ack store, `harness_acks.check`, `harness_review.py`, the push gate | `base_branch`, `forge` | master, github |
| The governed manifest, the self-protecting set, the state dir, `profile.toml` | `profile.rule_files` (additions) | as generic-harness §3.1 |
| New hooks: `AskUserQuestion`, `PermissionRequest`, `SubagentStart`; the recap header | none | — |

**Proposed generic rule text** (for `REQUIRED-READING`'s harness section;
rules only):

> **Unattended mode.** Only the owner enters it, with `tools/away.py
> <duration>` in a terminal; it ends at the time the flag states or on
> `away.py --back`. No agent writes the flag or runs either command. While
> it is on, a guarded act is refused rather than asked: do not retry or work
> around it; record it as an `ASK <OWNER>:` line (a subagent: in its
> handback) and continue with other work. A rule-file edit on a
> non-default branch outside the live checkout is allowed and logged; it
> binds nobody until the owner acknowledges it. In either mode, a push or
> merge whose ref changes a rule file relative to the base branch is refused
> until the owner has acknowledged that exact diff with
> `tools/harness_review.py`.

---

## 8. Where this slots into generic-harness §7

Step 1 has landed (#120), and step 2 (citation pinning) is in flight. The
feature splits into three PRs. The first only makes things stricter and
fixes today's stalls. The third is the only one that loosens anything, and
it lands after the gate that makes the loosening safe. Each PR adds hook
wiring to `.claude/settings.json`, so each needs Ola's fresh yes. LOC
figures are **estimates** of production lines as `CLAUDE.md` §2 counts them.

| # | Step | Placement | Files | Estimate |
|---|---|---|---|---|
| U0 | **Probes** (no production code). (a) A deny inside a background `@tester`: does the reason arrive and does the subagent carry on? This combines with R-B's `agent_type` probe. (b) Does `!cmd` at Ola's prompt get a TTY? (c) Do `.claude/current-task/` writes prompt in a session isolated in a worktree? | with step 8's probe, but moved earlier | a temporary logging hook (settings: fresh yes) | 0 |
| U1 | **Queue, don't wait.** `harness_mode.py` (flag, `decide`, queue), `away.py`, `guard_push` and `guard_governance` routed through `decide`, and the `AskUserQuestion` and `PermissionRequest` hooks, plus the recap header. Rule-file edits are **denied** in unattended mode in this step | **after step 3** (which rewrites the same hooks' docstrings), before step 4. It fixes a live cost now, and the manifest can stay in `guard_governance` until step 4 | `.claude/hooks/*`, `tools/away.py`, `tools/session_state.py`, `.claude/settings.json` | ~220 |
| U2 | **Acknowledgement gate.** `harness_acks.py`, `harness_review.py`, the gate in `guard_push` (both modes), the attended prompt-ack arm, and the recap's `REWORK` lines | **after step 4**, so that `governed()` and the gate share one manifest in `profile.toml` | `tools/harness_{acks,review}.py`, `.claude/hooks/guard_{push,governance}.py`, `.claude/settings.json` (PostToolUse `Edit\|Write`) | ~280 |
| U3 | **Authorship on the diff.** The unattended `pass` arm with the journal, the "not in force" context, and the `SubagentStart` hook | after U2; independent of steps 5-9 | `guard_governance.py`, a `SubagentStart` hook, `.claude/settings.json` | ~110 |
| — | The same code moves into the harness repo in step 6 and comes back in step 7 unchanged; the drift gate covers it | steps 6-7 | — | — |

Tests are extra and belong to `@tester`. They should cover a planted
unacknowledged rule diff, a stale base, a post-ack edit, a rejection, a
rule edit through a variable path, a deny inside a subagent, and an
unparseable flag.

---

## 9. Open questions for Ola

1. **Scope of the flag.** Should it be per repository (the git common dir,
   shared by that repository's worktrees) or per machine (`~/.claude/`)?
   *Recommend per repository.* Each project vendors its own harness
   version, and hooks should not read user-global state.
2. **How you enter the mode.** It would be only from your own terminal
   (`away.py 8h`), and the main session would ask you to run it when you say
   you are leaving, instead of starting `caffeinate` itself. *Recommend
   yes.* Your memory note on keep-awake would change to match; that note is
   yours to edit. U0(b) tells us whether `!away.py` at the prompt also works.
3. **Limits on duration.** The duration would be mandatory, capped at 72 h,
   with the same ×1.2 buffer as keep-awake. *Recommend yes.* An exhibition
   longer than 72 h means entering the mode again, remotely or before you
   leave.
4. **The self-protecting set.** Settings, hooks, the harness state,
   `away.py`, the review command and `profile.toml`, plus any rule file in
   the live checkout, would be denied at night even on a branch. *Recommend
   yes.* They change enforcement tonight, before your review can see them.
5. **Attended prompts count as acks.** Your yes at `guard_governance`'s
   prompt would record an ack of the resulting blob, so the review command
   is needed only for unattended edits and for edits the tripwire missed.
   *Recommend yes*, or daytime rule work would be reviewed twice. The
   stricter alternative is to review every rule diff before every push.
6. **Where acks are stored.** They would be untracked in the git common dir,
   rather than a tracked file or commit trailer that CI could check.
   *Recommend untracked for now.* A tracked record is as easy to forge on a
   shared machine and adds a governed file that agents would write around.
   Revisit it if you want a durable authorship trail for publication.
7. **History rewrites at night.** Should `amend` and `rebase` be denied even
   on commits that were never pushed? *Recommend yes for v1*: new commits
   are cheap, and a narrower rule means `guard_push` must know what is
   published.
8. **"Local" branch.** You said edits are allowed on local, non-default
   branches. Should that include a non-default branch that already has an
   open PR? *Recommend yes.* The push gate is what protects publication, and
   a follow-up PR round is the common case.

## Ola's rulings (2026-09-30)

All eight open questions are answered as recommended:

| # | Question | Ruling |
|---|---|---|
| 1 | Flag scope | Per repository |
| 2 | Who enters the mode | Ola only, from his own terminal (`tools/away.py <duration>`, which also starts `caffeinate`). When Ola says he is leaving, the main session asks him to run it instead of starting `caffeinate` itself |
| 3 | Duration | Mandatory, capped at 72 h, with the same ×1.2 buffer as `caffeinate` |
| 4 | Self-protecting set | Denied at night even on a branch |
| 5 | Daytime rule-file prompt | Ola's yes counts as the acknowledgement |
| 6 | Acknowledgement storage | Untracked for now; revisit if a durable authorship trail is wanted for publication |
| 7 | Amend and rebase at night | Denied, even on unpushed commits, in the first version |
| 8 | Branch with an open PR | Counts as a local branch; the push gate protects publication |
