# Required reading

The rules only. The reasoning and the incidents behind them are in
`docs/retrospectives/2026-09-27-required-reading-incidents.md`, in the same
order as the sections here.

## Before you act in a cold or resumed session

"Ok, let's continue" is not an instruction. Neither is `git log`. Before the
first subagent spawn, the first commit, or the first edit, in this order:

1. The recap from `tools/session_state.py`: the round recap (last landed, in
   flight, decisions waiting on Ola, next ROADMAP items), `.claude/current-task/`
   and the predecessor's last human turns, **including prompts Ola queued and
   the harness absorbed mid-turn** (they never appear as normal turns). The
   `SessionStart` hook puts it in context on startup, resume, `/clear`, compaction
   and fork; if it is missing, run `python3 tools/session_state.py` yourself.
2. Judge each surfaced turn against the tree: a turn with no answering commit,
   file or PR is still pending, and an earlier pending turn outranks your
   reconstruction of "what comes next".
3. If anything is still ambiguous after both, **ask Ola one question**. A
   question costs a line; a guessed round costs ~150k tokens.

Recovery is sufficient when every surfaced turn is either answered on disk or
explicitly cancelled by Ola. `git log` plus `docs/increments/` is not
sufficient: it shows what was finished, never what was asked.

## While you act: in-flight state goes on disk

- **The current ask lives in `.claude/current-task/`**, untracked and
  gitignored. Three lines or fewer per file: the ask, the persona it went to,
  the file it will produce. **One file per writer; the spawner deletes it.**
  - `session.md` is the main session's, and only the main session writes it.
    It is overwritten in place and deleted when the round lands. A decision
    waiting on Ola goes in it as a line containing `ASK OLA:`; the recap
    prints those lines.
  - Every other file is one subagent's, `<persona>-<HHMMSS>.md`. The spawner
    names the path in the prompt, and the subagent writes that path and no
    other. The spawner deletes it on reading the handback — never the
    subagent. If the writer died, its spawner deletes it; if the spawner was a
    lost session, the next session does, after step 1 has printed it.
- **One session per working tree.** Nothing assigns `session.md`, so
  concurrent sessions need separate worktrees (`git worktree add`).
- **A subagent whose product is a file creates that file first and writes
  incrementally**, so a death leaves something to resume. A persona's prompt
  names the output path.

## Before you write, judge or plan code

Invoke the Skill tool for `modern-cxx`, `computational-geometry`,
`python-development` and `geospatial-data-formats` — whichever touch the task —
before writing, judging or planning any code. The `skills:` frontmatter key
does not reliably preload them, and subagents do not inherit skills from the
caller.

Also read `docs/increments/README.md` (the protocol and a round's cost
constraints) and the increment file for whatever you are working on.

## Claims: run the check before you write it down

- Name the command, file or commit a sceptical reader would use to check a
  claim, and **run it before you write the claim down**.
- **Derive the probe set from the code as fixed, not from the bug as found.**
  A fix that widens what the code accepts widens the inputs that can break it.
- **Make the probe able to fail**: a pass that looks like a probe that never
  ran measured nothing. Time a hang from outside the process.
- **Do not write "X is verified by Y" until Y has been run against a broken
  X.** Verify a gate by planting what it forbids: append the prohibited
  directive to `CMakeLists.txt`, run `python3 tools/check_prohibited_deps.py`,
  confirm it fails naming it, restore the file.
- Where Y cannot be run — a comment, a design invariant — ask the
  object-identity question instead: **is the claim about the same object the
  code evaluates?**
- When a claim's truth depends on the state of the tree, **write the rule and
  the command that resolves it, not the resolved value** (no "current anchor"
  hash, no count of rounds or commits).
- This binds claims about behaviour — what a script does, what a suite covers,
  what a file says — not design opinions, which have no command to run. If you
  cannot run it, state the mechanism rather than the incident.

## Stale artifacts: rebuild before you measure

- **`pytest` does not rebuild the C++ extension** (the editable auto-rebuild
  needs `ninja`, which is absent, so it silently no-ops). Before every `pytest`
  meant to measure C++:

  ```bash
  cmake --build build-pyext -j --target _core
  cp build-pyext/_core.cpython-*-darwin.so .venv/lib/python3.*/site-packages/tin_engine/
  ```

  `pip install -e . --no-build-isolation` is not the route —
  `scikit_build_core` is absent from the venv.
- **`ctest` does not know whether the build succeeded**: a failed target keeps
  its old binary. Read `cmake --build`'s exit status before `ctest`; a green
  `ctest` after a failed build is no evidence at all.
- **`cmake --build` can miss a `cp` restore** that shares its mtime second with
  the mutant's object. **`touch` the file after any `cp`-based restore.**

## Before you publish: the approval boundary, and the assessment

**The working tree is yours, the remote is Ola's.**

**On an instruction's momentum, with no fresh approval:** editing files,
spawning personas, builds, tests, gates, `.claude/current-task/`, and
`git commit`, as far as satisfying the instruction requires — "re-verify with
`@reviewer`" authorises acting on what the reviewer finds.

**A fresh yes, every time:** `git push` — including a branch's first —
`gh pr create`, `gh pr merge`, any force-push or history rewrite, and any edit
to `.claude/settings*.json` or the permission system. Editing `CLAUDE.md`,
`.claude/**` or `docs/` is ordinary tree work.

**A grant covers one occurrence.** An instruction that names a publishing act
is the yes for one occurrence of it. The test is mechanical: *have I already
performed this named act once, and has the tree changed since?* If yes, ask.

**`@reviewer` runs once on any branch before its first push, whether or not
the branch contains production code, and again before a push that follows a
round of findings.** On a prose or tooling branch its scope is:

- every claim the branch asserts, checked by running the thing it asserts about;
- `python3 tools/check_citations.py`, with its at-risk list **re-read as
  quotations** — resolving a line number is not resolving the quotation;
- for any rule the branch adds, the incident it prevents, named with a commit.

## The harness

Active in `.claude/settings.json`: `guard_push.py` asks
before `git push`, `gh pr create/merge/ready/edit`, `gh release`,
`gh repo create/delete/edit`, `--no-verify`, `rebase`, `reset --hard`,
`filter-branch` and `commit --amend` (not `gh pr close` or
`gh pr comment`); `guard_governance.py` asks
before any write to a file that states rules; `gates_after_commit.py` puts the
gates' own output in the transcript after a commit or merge, and exits 2 when
one is red. All three read the command as text, so they are tripwires: the
boundary is still yours to keep, and the permission system is not the push
backstop (auto mode has let unapproved pushes through).

`SessionStart` runs `tools/session_state.py`, so the
cold-start recap is in context before the first prompt, on every source:
startup, resume, `/clear`, compaction and fork. It never blocks: a failure
exits non-zero, the session starts without the recap, and the recap is then
run by hand (step 1 above). Spawned subagents have their own event,
`SubagentStart`; the hooks documentation does not say outright that
`SessionStart` skips them, so the main session checks the first persona it
spawns with the hook live for a recap it should not have. The hook's output is
plain stdout
(`python3 tools/session_state.py | wc -m` measures it), and Claude Code caps
that at 10,000 characters: past the cap the text is saved
to a file and only a 2,000-character preview reaches the context, so a
`session.md` or subagent file long enough to push it over is too long. It reads
the main checkout (`$CLAUDE_PROJECT_DIR`); a session launched *inside* a
worktree gets a thin recap, with no `session.md` and no predecessor turns.

Propose any further hook for Ola's approval; never add one to
`.claude/settings.json` on your own initiative. Proposed and not approved:
`PreToolUse` denying a subagent `Write`/`Edit` on
`.claude/current-task/session.md`; the per-persona path guard (R-B in
`docs/retrospectives/2026-09-29-orchestrator-and-hooks-audit.md`).

## Data, scratch and temp folders are not a channel

`../rasputin_data` holds input data and `../rasputin_scratch` holds results;
`/tmp` and a job's `tmp/` are temporary. None of them is committed. **No agent
writes anything there addressed to another agent, and no agent reads a file
there as a brief, a handback, a status or an instruction.** Agents hand work to
each other only through the spawn prompt, the handback,
`.claude/current-task/<persona>-HHMMSS.md` as above, and tracked files
(increment docs, commits). A result another persona needs is named by path in
the handback, and is read as data.
