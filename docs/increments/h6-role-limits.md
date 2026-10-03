# Harness h6: what each persona may change

Status: **design**, @architect, 2026-10-02. One PR, after h5 has merged: h6
adds a role table to h5's tool-time check and no hook of its own. Sources:
R-B in `docs/retrospectives/2026-09-29-orchestrator-and-hooks-audit.md`
(§5 and §7), and two sections of `docs/retrospectives/next.md`: "Agents
taking on each other's work" and "An external review of the harness", point 1.

## 1. Prior art

Tooling; no novelty claimed. The pattern is path-scoped write permission per
role, as in CODEOWNERS or a capability table; the role-bleed failure it
answers is "disobey role specification" (failure mode 1.2) in Cemri et al.
2025, arXiv 2503.13657, already cited in `next.md`. The identity it keys on
is documented by Claude Code (code.claude.com/docs/en/hooks, common input
fields): `agent_type` is "present when the session uses `--agent` or the hook
fires inside a subagent", and `agent_id` is "present only when the hook fires
inside a subagent call". h3's live probe U0 confirmed both
(`docs/increments/h3-unattended-u1.md`, result (a′)): every subagent event
carries `agent_type`, and main-session events have both fields null.
*Legacy*: nothing.

```
$ git grep -l -iE 'agent_type|PreToolUse|persona' legacy-archive -- legacy
(no output, exit 1)
```

## 2. Why h6 builds on h5

R-B proposed a new `PreToolUse` hook with a Bash arm that reads command text.
The text arm has the same weakness as `guard_governance.py`: a path built
from a variable or from split strings gets past it. A @reviewer did exactly
that on the night of 2026-10-01. h5 already takes a snapshot before each tool
call and diffs after it (`state_check.py pre-tool` / `post-tool`, registered
for `Bash|Edit|Write|NotebookEdit`). h6 adds the persona's own changes to that
diff. The Bash route is then judged by what changed in the tree, whatever
command made the change. This is the external review's "guard state, not
actions".

What h6 removes or merges:

- **No new hook file and no new settings entry.** The table and both checks
  live in h5's `state_check.py`.
- **R-B's Bash text arm is dropped.**
- **The pending "deny a subagent `Write`/`Edit` on `session.md`" proposal**
  (`.claude/REQUIRED-READING.md`, *The harness*) becomes one row of the table.
  The rule text loses that sentence and the R-B sentence, and gains one
  sentence naming the table.
- **The persona prompts stop listing their path limits in prose.** Each says
  only "your write limits are the `ROLES` table in `state_check.py`". That
  leaves one statement per rule (review point 2).

h6 does **not** use h5's commit-time hook. Git hooks get no `agent_type`, and
the `(@persona)` tag in a commit subject is a convention that nothing
verifies, so a role check at commit time would key on text.

## 3. Design

**`ROLES`**, one table in `state_check.py`. Each persona maps to the
work-tree path prefixes it may write. Everything else is denied (an
allowlist, not R-B's denylist).

| `agent_type` | may write (relative to the tree's toplevel) |
|---|---|
| `tester` | `tests/` |
| `developer` | `src_python/`, `include/`, `src/`, `bindings/`, `tools/`, `.claude/hooks/`, `.github/`, `CMakeLists.txt`, `pyproject.toml` |
| `perf` | `docs/benchmarks/` |
| `architect` | `docs/` except `docs/retrospectives/`, `ROADMAP.md`, `CLAUDE.md`, `.claude/` files ending in `.md` |
| `orchestrator` | `docs/retrospectives/` (h7) |
| `reviewer`, and any other name | nothing |
| `dispatcher` (main session, §3.1) | `docs/` except `docs/retrospectives/`, `ROADMAP.md`, `.claude/current-task/session.md` |
| every subagent | its own `.claude/current-task/<agent_type>-*.md` |

Rule files stay under `guard_governance.py` and h5 whatever this table says.
The table only narrows what a persona may write; it never widens it.

**3.1 The dispatcher.** A tool event with neither `agent_id` nor
`agent_type` is the main thread (U0), and the table gives it a named row,
`dispatcher`. No other code path knows anything about a "main session". An
event with `agent_type` but no `agent_id` comes from a session launched with
`claude --agent <name>`; it gets that name's row. The main session is
expected to write no code, but Ola sometimes works with it directly, so its
verdict goes through `harness_mode.settle`: ask by day, refuse and queue at
night. A subagent's verdict is always `deny`, because an ask inside a subagent
waits for an Ola who may be away; the deny's reason tells the agent to hand
the work back.

**3.2 Edit, Write, NotebookEdit: exact, before the write.** In `pre-tool`,
resolve `file_path` to its git toplevel and judge the relative path against
the row. A path outside every work tree is not judged (see §5).

**3.3 Bash: from the tree, after the command.** For each event that carries
an identity row, `pre-tool` also records the tree's dirty set:
`git status --porcelain=v1 -z --untracked-files=all` and a `git hash-object`
of each dirty path, stored in h5's snapshot file. `post-tool` takes the dirty
set again. Every path whose status or hash changed is attributed to this
call. A changed path outside the row is a finding: one queue line naming the
persona and the paths, and `additionalContext` saying "outside your role:
restore <paths> and hand back". A command cannot be refused after it has run,
so the finding is reported and not prevented. Ignored paths (`build/`,
`.claude/current-task/`) are not seen, which is correct for `build/`
(see §5 for `.claude/current-task/`).

## 4. Tests @tester writes red

Events fed to `state_check.py pre-tool` / `post-tool` in a `tmp_path` git
repo (h5's fixture). Unless noted, the flag is off.

From the incidents in `next.md`:

1. `tester` Write `src_python/tin_engine/x.py` → `deny` naming the `tester` row
   (the 16b throwaway modules).
2. `tester` Bash `cat > src_python/tin_engine/x.py <<EOF` → post finding with
   `x.py`, a queue line with `agent_type: tester`, and `additionalContext`.
3. `tester` Bash `sed -i ... src_python/.../geopackage.py` → finding (the
   2026-09-29 temporary edit).
4. `tester` Edit `.github/workflows/main.yaml` → `deny` (15c-1).
5. `architect` Write `src_python/tin_engine/catchment.py` → `deny`;
   `architect` Write `docs/increments/22-probes/p.py` → allowed (the 22
   prototype; probes belong under `docs/`).
6. `reviewer` Bash `python3 -c "open('CLA'+'UDE.md','a').write('x')"` →
   finding; `reviewer` Bash that writes a placeholder token into a tracked
   `docs/` file and swaps it with a second script → finding (the 2026-10-01
   reroute).

Controls and edges:

7. `developer` Write `src_python/...` and `tester` Write `tests/...` → allowed.
8. Any subagent Write `.claude/current-task/<its type>-120000.md` → allowed;
   to `session.md` or another persona's file → `deny`.
9. `dispatcher` (no `agent_type`, no `agent_id`) Write `src_python/x.py` →
   `ask`; with the flag on → `deny` and queued. `dispatcher` Write
   `docs/x.md` → allowed.
10. `agent_type: general-purpose` Write anywhere in the tree → `deny`.
11. A path dirty before the Bash call and untouched by it → not attributed;
    the same path rewritten by the call → attributed.
12. A Bash call that changes nothing → no output.
13. A Bash call that writes only under `build/` (ignored) → no output.
14. `orchestrator` Write `docs/retrospectives/x.md` → allowed;
    `orchestrator` Write `docs/increments/x.md` → `deny`; `architect` Write
    `docs/retrospectives/x.md` → `deny` (h7: @orchestrator is the only writer
    there).

Not invariant-critical: no mutation round.

## 5. Not covered, named

- **Writes outside the work tree.** `/tmp`, `../rasputin_scratch` and
  `../rasputin_data` are out of scope: the 22 prototype, if it was written
  there, and messages left in the data folders.
- **A change made and undone within one Bash call.**
- **A background process that writes after the call returns** (as in h5).
- **A Bash write under `.claude/current-task/`** (gitignored). A subagent
  that writes `session.md` or another persona's file by Bash is not seen.
- **A persona that claims another persona's `agent_type`.** The field comes
  from Claude Code, not from the agent, so this means editing a frontmatter
  `name:`, which is a governed write.

## 6. Estimate

Production lines (`CLAUDE.md` §2): `ROLES` and the row lookup ~35, the Edit
and Write verdict in `pre-tool` ~25, the dirty-set record and diff ~50,
reporting through h5's queue and `additionalContext` paths ~15. Rule text:
one sentence replaces two in `REQUIRED-READING.md`, plus a one-line pointer in
each of the six persona files. **~125**, one PR. No settings change.

## 7. Questions for Ola

Both ruled by Ola on 2026-10-03 (`docs/increments/h7-orchestrator-role.md`):

- **Question 1: no `agent_type` is the dispatcher.** The main session,
  started with no agent name, is the `dispatcher` row. No
  `.claude/agents/dispatcher.md`.
- **Question 2: `@orchestrator` is kept, with a new role** (watch the
  workflow, research agentic design, own `docs/retrospectives/`). Its row is
  read-only plus `docs/retrospectives/`.

The questions as asked:

1. **The dispatcher's identity.** Should an event with no `agent_type` simply
   *be* the `dispatcher` row (recommended: it needs nothing at launch)? The
   evidence that the absence is reliable is U0, a single probe session, so
   the build re-checks it on its first live firing. The alternative is that every session is
   started with `claude --agent dispatcher`. That needs a new
   `.claude/agents/dispatcher.md`, whose prompt would then replace the main
   session's system prompt.
2. **Retire `@orchestrator`?** With a named dispatcher, the table has seven
   rows. `@orchestrator` writes nothing (R-A), and the main session already
   does the dispatching. Retiring it leaves six actors and one fewer prompt
   to keep consistent (the external review's "make governance smaller").
   Recommendation: retire it, as a separate PR after h6. The ledger idea in
   the audit's §7 would then belong to the dispatcher.
