# Harness h16: guard fixes, a line counter, a scratch copy, brief fixes

Status: Ola ruled on §7 on 2026-10-05 (all three defaults) and on the afternoon questions (last section). PR A: pushed as #185 (287 net production lines by `tools/count_loc.py bc01cd8 a6b966e`, against an estimate of 187). PR B, on `worktree-h16b`: red `9cf533d`, green `15c76f6` (78 net production lines by `tools/count_loc.py 98e31cd 15c76f6`, against an estimate of 60), PR A's head merged in as `b58ab57`, R1's *The harness* sentence written, `gh help` passed (red `9687c16`, green `084a6b3`); code review round 1 of PR B (round 6 below) asked for changes, and Ola chose option C: G3a kept for a single plain command, G3b dropped, G1's fetch rule and the merge guard widened; design `8554a3e`, red `94989dd`, green `a42d864` (77 net production lines by `tools/count_loc.py origin/master a42d864`). Code review round 2 of PR B (round 7 below) asked for changes; Ola ruled that PR B also closes the glued `gh api`/`curl` route, and then, on §7 question 4, the command-runner route past the push guard (last section). §2 G1 and G6 are amended and G7 is added for them, with §4, §6, §7 and R1 following; design `28b8292` and `1a1e5d3`, red `8c309a8`, green `fa9f3e1` (97 net production lines by `tools/count_loc.py origin/master fa9f3e1`). Code review round 3 of PR B (round 8 below) asked for changes; Ola ruled that PR B also closes the three routes it found (last section): §2 G7 is amended (every shell; a git or gh word under `parallel`) and G8 is added (a copy or move into a directory named without a trailing `/`), with §4, §6 and R1 following; `origin/master` merged in as `f86d3953`, red `2d44b462`, green `a0001d00` (109 net production lines by `tools/count_loc.py origin/master a0001d00`). Code review round 4 of PR B (round 9 below) asked for changes; Ola chose option A (last section): §2 G7 is amended (a runner in front of `parallel`, `watch` or `flock`, and the slow line of §7 question 5, which waits on Ola with default yes) and G4 is amended (a path judged also normalised), with §4, §6 and R1 following; writes to a whole governed directory go to §6. Ola answered §7 question 5 with its default on 2026-10-06 (fix the slow line); round 9's red step is `c83dfe0e` and its green step `80800a14`. Code review round 5 of PR B (round 10 below) asked for wording changes only: §2 G7's claim of linear time is corrected, and Ola ruled its suggested cap on runner and shell words goes to §6; the matching docstring fix is `85d39f5c`. Code review round 6 of PR B (round 11 below) approved it. Next, in order: a merge of `origin/master` by `@developer`, end to end (`git rev-list --count HEAD..origin/master` shows how far behind the branch is), then the push, on Ola's yes. `git log --oneline origin/master..HEAD` shows which of these have landed.

Ola approved the items on 2026-10-05 (the main session's summary of his
rulings, not his words). He said this is the last harness increment before
a freeze of a few days, so the design is kept to what each item needs and
§5 names what should drop out. What the labels mean:

| Label | What it is | PR |
|---|---|---|
| G1 | `guard_push.py` asks before `git fetch` into a named ref and before `git replace` writes | B |
| G2 | `guard_push.py` asks before a git or gh command it does not know (an alias from user config) | B |
| G3 | scratchpad repositories are ordinary (night proposal N3): G3a, file writes there; G3b (local git writes there) and G3c (the push half) dropped | B |
| G6 | the merge guard sees `gh` wherever `-R`/`--repo` stands, and `gh pr new`, `gh repo new` (review round 6); glued and clustered `gh api` and `curl` options (review round 7) | B |
| G7 | `guard_push.py` judges a git or gh command run by another program (`caffeinate git push`, `find … -exec git push`, `watch 'git push'`; §7 question 4), with a runner in front of `parallel`, `watch` or `flock` too (review round 9) | B |
| G4 | guard import shadowing: `tools/` last on `sys.path`, stdlib-named `tools/` files governed, a path judged also with `.` and `..` resolved (review round 9) | B |
| G5 | rule files written through the shell: one brief line | A |
| T1 | `tools/count_loc.py`, the one counter for `CLAUDE.md` §2 (night proposal N1) | A |
| T2 | `tools/scratch_copy.py <rev> <dir>` (night proposal N2), and evening P9 | A |
| T3 | `tools/brief.py`: note-file names that cannot collide; root Markdown files in `@architect`'s limit | A |
| R1 | rule lines: Monitor not `sleep`; a fallback kept for the last hour | A |

"The scratchpad" is the session's temporary directory,
`/private/tmp/claude-<uid>/<project>/<session>/scratchpad/`. "Governed"
means `guard_governance.py` asks before a write (`governed()` in that file).
Sources not on `master` are named by branch and commit:
`worktree-retro-1004c` holds the night retrospective
(`docs/retrospectives/2026-10-05-night.md`, a79f2d9) and the evening one
(`docs/retrospectives/2026-10-04-evening-merges-h11-h13-29.md`, c0e3451);
the brief said the evening one is on `master`, and it is not (`git branch -a
--contains c0e3451` lists only `worktree-retro-1004c`). `worktree-h12-design`
holds `docs/increments/h12-prose-fast-lane.md`.

## 1. Prior art: legacy and literature

Tooling; no novelty claimed.

*Literature.* The sources each item rests on, read 2026-10-05:

- git 2.55, `git-config(1)`, `alias.*`: "To avoid confusion and troubles
  with script usage, aliases that hide existing Git commands are ignored
  except for deprecated commands." So an alias can carry any name that is
  not a current command, and can take a deprecated one (checked here:
  `git -c alias.whatchanged='rev-parse --short HEAD' whatchanged` printed a
  hash). G2 rests on it.
- `git-config(1)`, `url.<base>.pushInsteadOf`: "Any URL that starts with
  this value will not be pushed to; instead, it will be rewritten to start
  with <base>, and the resulting URL will be pushed to." G3c rests on it.
- `git-replace(1)`: "Typing "git replace" without arguments, also lists all
  replace refs"; `-l`/`--list` lists. G1 treats those forms as reads.
- `python3(1)`, `-I`: "In isolated mode sys.path contains neither the
  script's directory nor the user's site-packages directory"; `-P`: "Don't
  automatically prepend a potentially unsafe path to sys.path such as [...]
  the script's directory". Python's *The initialization of the sys.path
  module search path*: "The first entry in the module search path is the
  directory that contains the input script". G4 rests on these.
- Anthropic, *Writing effective tools for agents* (cited by the night
  retrospective for N1 and N2) and *Making Claude Code more secure and
  autonomous with sandboxing* (for N3); not re-read here.

What differs: the h12 design scrubs user config
(`GIT_CONFIG_GLOBAL=/dev/null GIT_CONFIG_NOSYSTEM=1`) so that it reads only
the governed `.git/config`. That is right for a tool that must not be
misled (T1, T2, G3a's lookup) and wrong for a guard that must predict what
a command will do, because the command reads the user config the scrub
hides. G2 and G3c therefore do not scrub; §2 says where each applies.

*Legacy.* Nothing. `git grep -l -i -e 'insteadOf' -e 'sys.path.insert' -e
'git archive' -e 'count_loc' legacy-archive -- legacy` returned no file
(exit 1); the legacy tree had no harness.

## 2. The items

Each item: what changes, the incident it prevents, the red test, its size.
Sizes are production lines by `CLAUDE.md` §2 (net), estimated; tests are
excluded.

### G1. `git fetch` into a named ref, and `git replace`

**Change.** In `guard_push.py`'s `segment_why`:

- `git fetch` or `git pull` with a positional other than the first (the
  repository) that contains `:` asks: "fetch writes a named ref". That is a
  refspec with a destination (`. <sha>:refs/remotes/origin/master`,
  `+a:b`). `git fetch origin`, `git fetch -q origin master` pass, as today.
- `git replace` asks ("replace refs change what git reads for an object")
  unless it is a listing: no arguments, or `-l`/`--list` (with or without a
  pattern), or `--format=…` with those. Every other form (create, `-d`,
  `-f`, `--edit`, `--graft`, `--convert-graft-file`) asks.

**Incident.** h12 design review round 2 (recorded in dfe3ad7 on
`worktree-h12-design`): in a scratch repository,
`git fetch . "<sha>:refs/remotes/origin/master"` repointed `origin/master`
at a planted commit, and `git replace` made `git archive` serve a planted
file; neither was asked. h12's own text (§3.3, "These checks do not depend
on any guard change") says the guard change makes the routes harder, not
h12's checker correct; that is all G1 claims.

**Red test** (`tests/python/test_guard_push.py`, parametrised over argv):
asks for `git fetch . abc:refs/remotes/origin/master`,
`git fetch origin +master:refs/heads/x`, `git pull . a:b`,
`git replace HEAD HEAD~1`, `git replace -d x`, `git replace --graft a b`,
`git replace --edit a`; passes `git fetch origin`, `git fetch -q origin master`,
`git fetch --all`, `git replace`, `git replace -l`, `git replace --list 'a*'`.
Pinned false positive: `git fetch --depth 1 git@github.com:a/b.git` asks
(the URL is the second positional after `1`); `--depth=1` does not.

**Size.** About 12 lines.

**Green's departures** (`15c76f6`, `@developer`; both widen what asks, none
narrows it). The fetch check reads every word after the repository that does
not start with `-`, option values included, so a value can only move a
refspec later, never hide it. And `git_call` skips git's value-taking global
options (`-C`, `-c`, `--git-dir`, `--work-tree`, `--namespace`,
`--attr-source`) with their values when finding the subcommand, where it
skipped only `-C` and `-c` before.

**Amendment after review round 6 (Ola's option C).** A fetch can write a
named ref with no `src:dst` word on the line. Probed with git 2.55 in two
scratch repositories, each of these created a branch in the fetching one:
`git fetch --refm '+refs/heads/*:refs/heads/got/*' ../src master` (the
separate value is today taken for the repository, so the rule above reads
`../src master` and passes); `echo 'refs/heads/master:refs/heads/x' | git
fetch --st ../src`; `git -c remote.s.url=../src -c
'remote.s.fetch=+refs/heads/*:refs/heads/x/*' fetch s`, the same with the
key spelled `REMOTE.s.FETCH`; `git --config-env=remote.s.fetch=RF fetch s`
with `RF` holding the refspec; and `git -c … remote update s`. A
`remote.<x>.url` override alone re-points `refs/remotes/<x>/*` at planted
commits, which is G1's incident itself, and `url.<x>.insteadOf` does the
same through the URL. So, for a git command whose subcommand is `fetch` or
`pull`, or `remote` with first positional `update`, it also asks ("fetch
writes a named ref") when any of these holds:

- an argument whose part before `=` is at least `--ref` long and a prefix of
  `--refmap`, or at least `--st` long and a prefix of `--stdin`. Git takes
  any unambiguous prefix of a long option: with git 2.55, `--refm` and `--st`
  are the shortest it accepts for `fetch`, `--ref` for `pull` (`pull` has no
  `--stdin`; `--st` there is ambiguous and git refuses it). A prefix git
  refuses as ambiguous asks too, which costs nothing. `--refetch`, a
  different option, passes;
- among git's own options before the subcommand, a `-c <key>[=<value>]`,
  `--config-env <key>=<var>` or `--config-env=<key>=<var>` whose key,
  lower-cased, starts with `remote.` or `url.`. Git does not accept a glued
  `-c<key>` (`unknown option`), so that form needs no rule;
  `--config-env` joins `GIT_TAKES_ARG`, since git accepts its value as the
  next word;
- the line's text contains `GIT_CONFIG` (`GIT_CONFIG_COUNT`/`KEY_n`/`VALUE_n`,
  `GIT_CONFIG_PARAMETERS`, `GIT_CONFIG_GLOBAL`): `tools/shell_scan.py` strips
  leading assignments from the argv, so only the text shows them. A variable
  exported by an earlier command or a shell startup file is not seen (§6).

A `-c` with any other key (`git -c protocol.version=2 fetch origin`), and
`git remote update` with no override, pass as before. [Amended after review
round 7: the sentence that stood here, that such a key writes only what the
repository's own governed config says, was false; the next paragraph is the
correction.]

**Amendment after review round 7.** Probed with git 2.55 in two scratch
repositories, `src` holding a planted commit and `dst` fetching (setup and
results in this run's transcript, no commit). Each of these wrote a named
ref in `dst` (`refs/heads/got/master`), with `F` a file holding
`[remote "s"] url = …/src` and `fetch = +refs/heads/*:refs/heads/got/*`:
`git -c include.path=F fetch s`, the same as `-c INCLUDE.PATH=F`,
`git -c includeIf.onbranch:master.path=F fetch s`,
`git -c includeIf.gitdir:<dir>/.path=F fetch s`,
`git --config-env=include.path=V fetch s` (`V` holding F's path),
`HOME=<dir> git fetch s` (F as `<dir>/.gitconfig`),
`XDG_CONFIG_HOME=<dir> git fetch s` (F as `<dir>/git/config`), and
`GIT_CONFIG_GLOBAL=F` or `GIT_CONFIG_SYSTEM=F` before `git fetch s`. And,
with the repository's `origin` at its real ssh URL
(`git@github.com:expertanalytics/rasputin.git`),
`git -c core.sshCommand=<program> fetch origin`, the program running
`git-upload-pack <dir>/src`, re-pointed `refs/remotes/origin/master`,
`…/other` and `…/HEAD` at the planted commit: G1's incident, through the
transport. `GIT_SSH=<program> git fetch origin` exits 0 the same way.
`git -c fetch.bundleURI=file://<bundle> fetch <dir>/src master` wrote
`refs/bundles/heads/master`.

The probe set comes from `git help --config` (git 2.55), every key that
bears on where a fetch goes or what it writes, not only the reported route:
`remote.*`, `url.*`, `include.path`, `includeIf.<condition>.path`,
`branch.<name>.*`, `fetch.*`, `bundle.*`, `transfer.*`, `core.sshCommand`,
`core.gitProxy`, `core.askPass`, `core.alternateRefsCommand`, `ssh.variant`,
`protocol.*`, `http.proxy*`, `credential.*`, `uploadpack.*`. Their fate:

- **Ask** (added to the key list above, compared lower-cased): a key
  starting `include.` or `includeif.` (a whole file of config, so any of
  the rest); exactly `core.sshcommand` (the transport for `origin`, an ssh
  URL, so it chooses what `refs/remotes/origin/*` receive); exactly
  `fetch.bundleuri` (a URL whose bundle is unpacked into `refs/bundles/*`).
  The text check widens from `GIT_CONFIG` to the regular expression
  `GIT_CONFIG|GIT_SSH|\b(HOME|XDG_CONFIG_HOME)=` on the line's text:
  `GIT_SSH` covers `GIT_SSH_COMMAND`; `\bHOME=` does not match
  `JAVA_HOME=`, since `_` is a word character.
- **Pass, pinned:** `branch.<name>.remote` and `.merge` choose which remote
  and which ref a plain `git pull` takes, which the command line already
  does unasked (`git pull . other`); a URL there fetches into `FETCH_HEAD`
  only. Other `fetch.*` keys (`fetch.prune`, `fetch.pruneTags`) delete stale
  tracking refs as `--prune` does on the command line, which passes; they
  never point a ref at a new commit. `core.gitProxy` serves only `git://`
  URLs, and `http.proxy*` and `credential.*` only `http(s)://`; `origin` is
  ssh, so they need a `remote.`/`url.` override, which asks.
  `protocol.*`, `ssh.variant` and `transfer.*` change how, not where.
  `uploadpack.*` and `bundle.*` act only with a local path remote or a
  bundle URI, each of which needs a key that asks or writes only
  `FETCH_HEAD`.
- **Not closed, §6:** `core.askPass`, `core.alternateRefsCommand`,
  `uploadpack.packObjectsHook`, `core.fsmonitor`, `core.hooksPath` and
  `credential.helper` run a program, on commands other than fetch too, so
  asking on fetch would not close them.

Two persistent routes, found while probing, that are not one-command
config: git 2.55 still reads the legacy remote files. With
`<git-dir>/remotes/s` holding `URL: …/src` and
`Pull: refs/heads/master:refs/heads/got/viaremotes`, a plain `git fetch s`
wrote `refs/heads/got/viaremotes`; with `<git-dir>/branches/s` holding
`…/src#master`, it wrote `refs/heads/s`. Neither path is governed today
(`governed('.git/remotes/s')` and `governed('.git/branches/s')` return
False, checked on `a42d864`). So `GOVERNED_PREFIXES` in
`guard_governance.py` gains `.git/remotes/` and `.git/branches/`, as
`.git/hooks/` is there. The user's own `~/.gitconfig` and
`~/.config/git/config` are not governed either (§6).

**Red test, round 7 amendment** (`test_guard_push.py`, both modes, reason
"fetch writes a named ref"): asks for `git -c include.path=/x/f fetch s`,
`git -c INCLUDE.PATH=/x/f fetch s`, `git -c includeIf.onbranch:master.path=/x/f
fetch s`, `git -c includeif.gitdir:/x/.path=/x/f pull s`,
`git --config-env=include.path=F fetch s`, `git --config-env include.path=F
remote update s`, `git -c core.sshCommand=/x/p fetch origin`,
`git -c fetch.bundleURI=file:///x/b fetch origin`, `HOME=/x git fetch
origin`, `XDG_CONFIG_HOME=/x git fetch origin`, `env HOME=/x git fetch
origin`, `GIT_SSH_COMMAND=/x/p git fetch origin`, `GIT_SSH=/x/p git pull
origin`; passes `git -c fetch.prune=true fetch origin`,
`git -c branch.master.remote=origin pull`, `JAVA_HOME=/x git fetch origin`,
`git -c include.path=/x/f status` (not a fetch), `git -c
core.sshCommand=/x/p log -1`. `test_guard_governance.py`: `governed()` is
True for `.git/remotes/s`, `.git/branches/s`,
`/abs/repo/.git/remotes/origin`; False for `docs/remotes/s`. About 3
lines.

**Red test, amendment** (`test_guard_push.py`, both modes, reason "fetch
writes a named ref"): asks for `git fetch --refmap=+a:refs/heads/x o master`,
`git fetch --refmap '+a:refs/heads/x' o master`, `git fetch --refm=+a:b o
master`, `git pull --ref=+a:b o master`, `git fetch --stdin o`, `git fetch
--st o`, `git -c remote.s.fetch=+a:b fetch s`, `git -c REMOTE.s.FETCH=+a:b
fetch s`, `git -c remote.s.url=/x fetch s`, `git -c url./x.insteadOf=https://github.com/
fetch origin`, `git --config-env=remote.s.fetch=RF fetch s`, `git
--config-env remote.s.fetch=RF fetch s`, `git -c remote.s.fetch=a:b pull s`,
`git -c remote.s.fetch=a:b remote update s`, and `GIT_CONFIG_COUNT=1
GIT_CONFIG_KEY_0=remote.s.fetch GIT_CONFIG_VALUE_0=a:b git fetch s`; passes
`git fetch --refetch origin`, `git -c protocol.version=2 fetch origin`,
`git remote update`, `git -c remote.s.fetch=a:b status`. About 8 lines.

### G2. Git and gh aliases: user config changes what the guard sees

**Change.** `guard_push.py` judges a git command by its subcommand word, so
an alias hides what runs. Probe, run on this branch's base (bc01cd8):
`publishes(['git', '-c', 'alias.p=push', 'p', 'origin'])` and
`publishes(['git', 'p', 'origin'])` both return `[]`, and `segment_why`
returns `None`: an alias written into `~/.gitconfig` (not governed) or given
with `-c` pushes unasked. The fix: a git subcommand that is not a current
git command asks, with the reason "a git alias or a command the guard does
not know; it cannot see what it runs". "Current" is the output of
`git --list-cmds=main` less `git --list-cmds=deprecated` (181 and 2 names
with git 2.55), read once per hook run. Those two calls are not scrubbed:
they list commands, which no config changes. The same for gh: `words[1]`
outside a fixed set of gh's top-level commands asks. The set is what
`gh help` lists with gh 2.101 less its alias `co` (which a user can
redefine): `auth browse codespace discussion gist issue org pr project
release repo skill cache run workflow agent-task alias api attestation
completion config copilot extension gpg-key label licenses preview ruleset
search secret ssh-key status variable`, plus `help`, which only reads
(Ola's ruling of 2026-10-05, quoted in `9687c16`'s test comment; green
`084a6b3`). A gh extension is then asked about too. Since review round 6
the gh word judged is the first one after gh's options (G6 below), so
`gh -R o/r pm 12` asks as `gh pm 12` does.

**Incident.** h12 design review round 3 (recorded in 2d77b6c on
`worktree-h12-design`): a user-level config file changed what git resolved
(`url.<x>.insteadOf` sent `ls-remote` and `fetch` elsewhere) with no
guarded command. That round's subject was the h12 checker, not the guards;
the alias route above is the same class turned on `guard_push.py`, found
while writing this design, with no incident of its own.

**Red test** (`test_guard_push.py`): asks for `git -c alias.p=push p origin`,
`git p origin`, `git whatchanged` (deprecated, so aliasable), `gh pm 12`,
`gh co 12`; passes `git status`, `git log -1`, `git worktree list`,
`gh pr view 12`, `gh api repos/x` (a GET). A test that `git --list-cmds=main`
contains `push`, so a git without `--list-cmds` fails loudly rather than
asking on everything.

**Size.** About 14 lines.

### G3. Scratchpad repositories are ordinary (N3)

Night retrospective §6, N3 (a79f2d9). Three parts were designed; after
review round 6 Ola chose option C (last section), so only G3a ships, and
narrowed.

**G3a, governance.** A path is exempt when it is absolute and its real path
(`os.path.realpath`, so `/tmp` resolves to `/private/tmp` and a symlink to
the repository's `.git` is followed) lies under a scratchpad: the pattern
`^/private/tmp/claude-\d+/[^/]+/[^/]+/scratchpad(/|$)`. Any session's
scratchpad, not only the current one: all are temporary, and nothing the
harness reads lives there. A relative path is judged as today, since the
guard does not track `cd`; so is a path built from an expansion
(`$SP/x/CLAUDE.md`), because the guard judges its static part
(`/x/CLAUDE.md`), which is not under a scratchpad.

Where the exemption applies (amended after review round 6):

- **Edit, Write, NotebookEdit:** always, as before. One tool call writes one
  path, and a link made by any earlier command is already followed by
  `realpath`.
- **Bash:** only when the line is a single plain command:
  `shell_scan.parse(command)` returns exactly one `Simple`, and that
  `Simple` has `program is None`, `script is None` and `unknown is False`.
  So a second command joined by `&&`, `||`, `;`, `|`, `&` or a newline, a
  `$(…)`, backticks or `<(…)` (each parses to a command of its own),
  `sh -c '…'` (its text parses to a second command), and an interpreter
  running a program or a script each void the exemption, and every target
  of that line is judged as if it did not exist. Checked with
  `shell_scan.parse` on this branch: `ln -sfn … && echo … > …`,
  `echo $(true) > …`, `sh -c 'echo > …'` and `printf x | tee …` give two
  commands; `python3 -c "open('…', 'w')"` gives one, with `program` set. A
  lone command in parentheses or followed by `&` parses to one command and
  keeps the exemption, which is harmless: it is still one command. A
  heredoc body is text, so `cat >> <pad>/r/.git/config <<'EOF'` is still one
  plain command. A line the parser cannot read is judged by
  its text, as before (`text_hits`, which never consulted the exemption).
  Why: the guard judges targets before the line runs, so
  `ln -sfn <real repo> <pad>/x && echo … > <pad>/x/CLAUDE.md` resolves
  `<pad>/x` while it does not yet exist, and passed (review round 6,
  finding 2). One plain command cannot make a link and then write through
  it, except a program that does both itself, which is the residual below.

Interface: `def governed(path: str, *, scratch_exempt: bool = True) ->
bool`. The scratchpad check runs first, only when `scratch_exempt` is true;
the rules below it are unchanged. The Edit/Write arm calls `governed(path)`;
`judge_bash` computes `plain` once per line and calls
`governed(shell_scan.static(t), scratch_exempt=plain)`. `tools/scratchpad.py`
stays (only `guard_governance.py` imports it now), and so does its
`GOVERNED` entry; its docstring drops "G3b".

Residual (also §6): a single command that itself makes a link and writes
through it (an archive extracted with a link and a file under it, a
compiled program), and a background job started earlier that re-points a
link between the hook's check and the write.

**G3b, local git writes in a scratch repository: dropped** (option C). It
passed `config`, `remote`, `symbolic-ref`, `update-ref`, fetch and replace
writes given as `git -C <absolute scratch dir>`. Review round 6 found that
it matched only exact spellings of the options that send git elsewhere
(`git -C <pad>/r config --glo …` wrote the global config;
`-f/path/to/real/.git/config`, glued, wrote the real repository's), and
that a line could swap the scratch repository's `.git` for a link before
the write (finding 2). What goes: `in_scratch`, `SCRATCH_VOIDS`,
`guard_push.py`'s import of `scratchpad`, `os` and `pairwise` if nothing
else uses them, `segment_why`'s exemption branch, and the G3b sentences in
its module docstring. `segment_why` keeps its `text` parameter for G1's
`GIT_CONFIG` check. A local git write in a scratch repository now asks like
any other; the routes that remain without a question are a file write to
`<pad>/r/.git/config` (G3a) and `git -c key=value` for a single command.

**G3c, a push to a scratch repository: dropped** (Ola's ruling 3, first
rulings section). A push always asks.

**Incident.** The night of 2026-10-04/05, §2 of the night retrospective
(a79f2d9): four of seven `guard_unattended` refusals were scratch
repositories in the scratchpad (22:58:22 `git config user.*`,
`git remote add` and a push to a local bare repository, agent `a5b4f008`;
23:04:28 an append of `[url] insteadOf` to a scratch `.git/config` and
23:05:27 a push to a `file://` scratch remote, agent `a3cc5c9a`), and two of
them were then reached by another route. A refusal leaves no commit; the
retrospective's table is the record.

**Red test, as first shipped** (`9cf533d`). `test_guard_governance.py`: a
Write of `<scratchpad>/r/.git/config` passes; of
`<scratchpad>/copy/CLAUDE.md` passes; of `<scratchpad>/link/.git/config`
where `link/.git` is a symlink to a real repository's `.git` asks; of
`.git/config` (relative) asks. A Bash `cat >> <scratchpad>/r/.git/config`
passes. These all stay, unchanged. `test_guard_push.py` had G3b's block,
whose fate is below.

**Tests after option C.**

*Stay, unchanged:* the G3a tests above (`test_a_file_under_a_scratchpad_is_not_governed`,
`test_a_shell_write_under_a_scratchpad_is_not_governed`,
`test_a_git_config_that_resolves_outside_the_scratchpad_asks`),
`test_scratchpad.py`, and `harness_fixtures.py`'s `point_scratchpad` and
its `tools/scratchpad.py` entry in `COPIED`.

*Go* (`test_guard_push.py`): `G3B_PASSED` and
`test_a_local_write_in_a_scratch_repository_is_silent`. Their rows do not
vanish: they move into the asked table below, where they are red until G3b
is removed.

*Stay, in a changed role:* `G3B_ASKED`,
`test_a_scratch_write_that_reaches_outside_still_asks`, the `scratch`
fixture and `make_plain_repo`. Renamed for what they now pin (for example
`G3_ASKED`, `test_a_git_write_in_a_scratch_repository_asks`), with the
comment block saying G3b was dropped by option C, and gaining:

- every former `G3B_PASSED` row with its reason (`CONFIG`, `REMOTE`,
  "update-ref moves a ref directly", "symbolic-ref rewrites a symbolic
  ref", "fetch writes a named ref", "replace refs change what git reads for
  an object"), except `git config --file {inside}/.git/config user.name x`
  and `cd {pad} && git -C {inside} config user.name x`, which ask too and
  are kept as rows;
- review round 6's spellings: `git -C {inside} config --glo user.name x`
  and `git -C {inside} config -f{outside}/.git/config user.name x` (both
  `CONFIG`).

Run in both modes `off` and `on`, as the rest of the file does (`on`
expects a queued `deny` carrying the same reason).

*New* (`test_guard_governance.py`, Bash events, both modes, each expecting
an `ask` naming the target; `{pad}` the fixture's scratchpad, `{repo}` the
copy):

- `ln -sfn {repo} {pad}/x && echo hi > {pad}/x/CLAUDE.md` (with `{pad}/x`
  absent before the hook runs: finding 2's route);
- `echo hi > {pad}/copy/CLAUDE.md; true`;
- `printf x | tee {pad}/copy/CLAUDE.md`;
- `echo $(true) > {pad}/copy/CLAUDE.md`;
- `sh -c 'echo hi > {pad}/copy/CLAUDE.md'`;
- `python3 -c "open('{pad}/copy/CLAUDE.md', 'w')"`.

And silent, both modes: `echo hi > {pad}/copy/CLAUDE.md`,
`cp {repo}/notes.txt {pad}/copy/CLAUDE.md`, and the existing heredoc
append. A unit test of `governed(p, scratch_exempt=False)` on a scratchpad
`CLAUDE.md` returning True pins the keyword.

**Size.** G3a about 9 lines (the exemption, and the plain-command test in
`judge_bash`); G3b's removal about −26 against PR B's green.

### G6. The merge guard and gh's `-R`/`--repo` (review round 6, finding 4)

**Change.** `publishes()` reads gh's group and verb as `words[1:3]`, and
`segment_why`'s G2 check skips any `words[1]` starting with `-`. With
gh 2.101 each of `gh -R a/b pr merge --help`, `gh pr -R a/b merge --help`,
`gh --repo=a/b pr merge --help` and `gh -Ra/b pr merge --help` printed
`gh pr merge`'s help, so each spelling reaches `pr merge`, and none asks
today. Also pre-existing: `gh pr new` is an alias of `gh pr create`, and
`gh repo new` of `gh repo create` (their `--help` lists them under
ALIASES); neither asks today. (`gh release new` asks, as all of `gh
release` does.)

One helper, used by `publishes`, the G2 check and the `gh api` check:

    def gh_words(words: list[str]) -> list[str] | None
        # None unless words[0]'s basename is gh; else words[1:] less every
        # word starting with "-", and less the word after an exact "-R" or
        # "--repo" (their glued forms "-R<x>" and "--repo=<x>" are one word).

The group is its first word and the verb its second. The PR writes become
`create`, `new`, `merge`, `ready`, `edit`, `update-branch`; the repo writes
`create`, `new`, `delete`, `edit`. On the groups the guard judges, `-R`/`--repo` is the only flag that takes a
value (`gh pr --help`, `gh release --help`, `gh repo --help`, and `gh --help`
for the top level, with gh 2.101), so no other value can be taken for the
verb. The text rules for a line the parser cannot read widen the same
way: `\bgh\b[^|;&]*\bpr\b[^|;&]*\b(create|new|merge|ready|edit|update-branch)\b`
and `\bgh\b[^|;&]*\b(release|repo\b[^|;&]*\b(create|new|delete|edit))\b`.
These can ask about an unreadable line that only names those words
(`gh pr list --search merge "`); that is the text rules' usual trade.

**Red test** (`test_guard_push.py`, both modes): asks, with "gh pr changes
a pull request", for `gh -R o/r pr merge 12`, `gh --repo o/r pr merge 12`,
`gh --repo=o/r pr merge 12`, `gh -Ro/r pr merge 12`, `gh pr -R o/r merge
12`, `gh pr new`, `gh -R o/r pr new`, and the unreadable `gh -R o/r pr
merge 12 "`; with "gh publishes or alters the repo" for `gh repo new x` and
`gh -R o/r release create v1`; with the G2 reason for `gh -R o/r pm 12`.
Silent: `gh -R o/r pr view 12`, `gh pr -R o/r view 12`, `gh --version`.

**Size.** About 8 lines.

**Amendment after review round 7 (Ola's ruling, last section): glued
options to `gh api` and `curl`.** `segment_why`'s `gh api` check (since h3)
counts a field only for an exact `-f`, `-F`, `--field`, `--raw-field` or
`--input` (or one of those before `=`), and its `curl` check counts data
only for an exact `-d`, `-F`, `--form`, `--json`, `-T`, `--upload-file` or a
word starting `--data`. Glued and clustered short options get past both, so
`gh api graphql -fquery='mutation { mergePullRequest(…) }'` merges a pull
request unasked, past the merge guard. On `a42d864`, `segment_why` returns
None for `gh api -fquery=x graphql`, `gh api -iXPUT
repos/o/r/pulls/1/merge`, `curl -sd x https://api.github.com/x`,
`curl --form-string a=b https://api.github.com/x` and `curl -sXPUT
https://api.github.com/repos/o/r/pulls/1/merge` (a `PUT` to that endpoint
merges with no body).

What the tools accept, probed:

- gh 2.101 (`gh api … graphql --help`, which exits 1 on a flag it cannot
  parse): `-ifquery=x`, `-iFquery=x`, `-iXPUT`, `-iiXPUT`, `-Xput`, `-X=PUT`
  parse; `--fie`, `--meth` (abbreviated long options) and `-iz` are
  refused. `gh api --help` lists the short flags `-F -H -X -f -i -p -q -t`,
  of which only `-i` (`--include`) takes no value, so a cluster is some `i`s
  followed by one value-taking flag and its glued value.
- curl 8.7.1, against a local listener: `-sd x`, `-d@f` and `-sFa=b` sent a
  `POST`; `-Tf`, `-sTf` and `-sXPUT` a `PUT`; `--form-string a=b` and
  `--json '{}'` a `POST`. `--requ`, `--data-b` and `--upl` are refused
  (curl takes no abbreviated long option). `--expand-data x` parses (curl
  has an `--expand-` form of each option).

**Change.** For `gh api`, over the words after the first `api` (as today):

- a word that starts with one `-` loses the run of `i`s right after it when
  more follows (`-iXPUT` → `-XPUT`, `-ifq=x` → `-fq=x`; `-i` and `-ii`
  stay), before both checks;
- a field is, besides today's forms, any such word starting `-f` or `-F`
  (`-fquery=x`, `-Fquery=@f`).

`method()` already reads `-XPUT` and `-X PUT`, so the normalised words need
nothing more. `-X=PUT` reads as `=PUT`, not `GET`, so it asks, as it
should; `-X=GET` with a field then asks too, a pinned false positive.

For `curl`, over all its words (as today), a word is data when:

- it starts `--` and, after a leading `--expand-` is cut to `--`, starts
  with `--data`, `--form` (so `--form-string`, and `--form-escape`, a
  harmless false positive), `--json` or `--upload-file`;
- or it is a short-option cluster, a word starting with one `-` and at
  least two characters long, that contains `d`, `F` or `T` anywhere.

And the method is also read from a cluster containing `X`: the text after
its first `X`, or the next word when `X` ends it (`-sXPUT`, `-sX PUT`), and
from `--expand-request`. "Anywhere" over-asks on a glued value that holds
one of those letters; pinned false positives: `curl -o/tmp/data.json
https://github.com/x` (`d`), `curl -HContent-Type:x https://github.com/x`
(`T`). These must still pass,
both modes: `curl -fsSL https://github.com/x`, `curl -sI
https://github.com/x`, `curl -s -o out https://api.github.com/x`,
`curl -sXGET https://api.github.com/x`, `curl -sX GET
https://api.github.com/x`; and `gh api repos/x`, `gh api -i repos/x`,
`gh api --paginate repos/o/r/pulls`, `gh api -q .name repos/x`,
`gh api -iXGET repos/x`, `gh api -X GET search/issues -f q=x` (as today).

**Red test, round 7 amendment** (`test_guard_push.py`, both modes): asks,
with "gh api with a writing method changes the forge", for `gh api graphql
-fquery=x`, `gh api graphql -Fquery=@f`, `gh api -ifquery=x graphql`,
`gh api -iFquery=x graphql`, `gh api -iXPUT repos/o/r/pulls/1/merge`,
`gh api -iiXPUT repos/o/r/pulls/1/merge`, `gh -R o/r api -fquery=x
graphql`; with "curl with a writing method to the forge" for `curl -d@f
https://api.github.com/graphql`, `curl -sd x https://api.github.com/x`,
`curl -Tf https://api.github.com/x`, `curl -sTf https://api.github.com/x`,
`curl -sFa=b https://api.github.com/x`, `curl --form-string a=b
https://api.github.com/x`, `curl --expand-data x https://api.github.com/x`,
`curl -sXPUT https://api.github.com/repos/o/r/pulls/1/merge`,
`curl -sX PUT https://api.github.com/repos/o/r/pulls/1/merge`,
`curl --expand-request PUT https://api.github.com/x`. Silent: the pass list
above. About 8 lines.

### G7. Commands run by another program (§7 question 4, Ola's ruling)

**Change.** `tools/shell_scan.py` strips ten wrappers (`env`, `nohup`,
`time`, `timeout`, `nice`, `command`, `exec`, `xargs`, `script`, `sudo`);
any other program that runs its arguments as a command hides that command
from `guard_push.py`. A helper in `guard_push.py`,

    def runs(words: list[str]) -> list[tuple[list[str], bool]]
        # (argv, bare) pairs: the command itself (bare False), and the
        # commands its later words run, by the three rules below.

and `main` judges every pair with `publishes` and `segment_why`, on both
paths (parsed simples, and `tokens(segment)` for a line the parser cannot
read). For each word after the first:

- **(a) A bare `git` or `gh` word**, by basename (`/usr/bin/git` too): the
  words from it to the end, with `bare` True. Any program, no list: a list
  of runners is never complete (`caffeinate`, `stdbuf`, `arch`, `flock`,
  `find -exec`, `parallel`, `watch`, and `taskset`, `unbuffer`, `ionice`
  beyond them), and a later `git` or `gh` word is what they share. A bare
  tail never gives G2's unknown-command reason: `grep git file` would
  otherwise read as `git file` and ask. Every other reason applies.
  [Amended after review round 8.] **Except under `parallel`**: when the
  command's own program is `parallel`, the tail is not bare, so it keeps
  G2's reason. `parallel` builds each command from its template and its
  inputs (`parallel git ::: push` runs `git push`; `parallel ::: git :::
  push` and `parallel git {} ::: push` do too), so the tail as written
  (`git ::: push`) is not the command that runs, and the guard cannot see
  it. A false positive is an ask on a harmless line such as `parallel git
  ::: status` or `parallel grep git ::: a` (pinned below); expanding the
  inputs exactly (GNU parallel's `{}` strings, `:::+`, `::::` input files)
  would cost more lines than the route is worth.
- **(b) A shell word**, basename in `shell_scan.SHELLS`: the words from it to
  the end, joined with `shlex.join` and parsed by `shell_scan.parse`, which
  already reads a shell's `-c` program as commands; each command it returns
  is judged in full (G2 included) and itself goes through `runs`. So
  `caffeinate sh -c 'cd x && git push'` asks. [Amended after review round
  8: the round found the list was `sh`, `bash` and `zsh` only, copied into
  `guard_push.py`.] The list is one set, `shell_scan.SHELLS`, which
  `guard_push.py` reads rather than copies, and `shell_scan.PROGRAM_FLAG`
  maps every name in it to `c`. It is derived from the host, not guessed:
  `/etc/shells` and `ls /bin/*sh` on Ola's Mac both give `bash`, `csh`,
  `dash`, `ksh`, `sh`, `tcsh`, `zsh` (run for this design); to those it
  adds the common shells that take `-c` and are not installed there:
  `rbash` (in Ubuntu's `/etc/shells`), `fish`, `mksh`, `ash`, `yash`. Each takes its program with `-c`, and the
  parser reads a simple command in that program whatever the shell's
  grammar. One set serves both guards: the first word (`dash -c 'git
  push'`, read by `shell_scan.parse` as a nested line) and a later word
  (`caffeinate dash -c 'git push'`, rule b); the governance guard reads the
  same nested line, so `dash -c 'cp x CLAUDE.md'` asks too.
- **(c) A quoted command line**, only when the command's own program
  (basename of the first word) is `watch`, `parallel` or `flock`: each later
  word containing whitespace is parsed by `shell_scan.parse`, and each
  command it returns is judged in full and goes through `runs`. These three
  take a command as one string (`watch 'git push'`, `parallel 'git push'
  ::: a`, `flock /tmp/l -c 'git push'`). Not every quoted word: then
  `grep -rn "git push" .claude/` and `git commit -m "git push is guarded"`
  would ask, both common in harness work and each a refusal while
  unattended.

A word that `shell_scan.parse` cannot read adds nothing (the outer line is
judged as before); with `shell_scan` missing, (b) and (c) add nothing. The
recursion ends because each parse is of a strictly shorter text.

The runner set, derived from the probes (§6, added after review round 7;
each passed through the hook on `a42d864` and on this branch's head
`28b8292`): as words before `git` (rule a) `find … -exec`, `caffeinate`,
`stdbuf`, `watch`, `flock`, `parallel`, `arch`; as a quoted string (rule c)
`watch`, `parallel`, `flock -c`; through a shell (rule b) `caffeinate sh -c`
and `find … -exec sh -c`. Probed with a scratch prototype of the three rules
(about 16 net lines, in this run's scratchpad, removed): every must-ask row
below asked with the reason given, every must-pass row passed, and the
pinned rows behaved as pinned; with the hook as on `28b8292`, all 22
must-ask rows passed.

**Incident.** Review round 7 (later item: `find -exec`), widened in this
design's §6 to the seven runners; no push was made that way.

**Red test** (`test_guard_push.py`, both modes). Asks, with the reason in
brackets:

- [PUSH] `find . -maxdepth 0 -exec git push origin HEAD \;`,
  `caffeinate -i git push`, `stdbuf -o0 git push`, `watch -n1 git push`,
  `flock /tmp/l git push`, `parallel git push ::: a`, `arch -arm64 git push`,
  `caffeinate -i /usr/bin/git push`, `nohup caffeinate stdbuf -o0 git push`,
  `watch -n1 'git push'`, `parallel 'git push' ::: a`,
  `flock /tmp/l -c 'git push'`, `caffeinate sh -c 'git push'`,
  `caffeinate sh -c 'cd x && git push'`,
  `find . -maxdepth 0 -exec sh -c 'git push' \;`,
  `watch 'caffeinate git push'`;
- [PR] `caffeinate -i gh pr merge 12`, `caffeinate -i gh -R o/r pr merge 12`;
- ["gh api with a writing method changes the forge"]
  `caffeinate -i gh api -X PUT repos/o/r/pulls/1/merge`;
- [FETCH] `caffeinate -i git -c remote.s.url=/x fetch s`;
- ["update-ref moves a ref directly"]
  `caffeinate -i git update-ref refs/heads/x HEAD`;
- [UNKNOWN, from rule b] `caffeinate sh -c 'git p origin'`.

Passes (false-positive controls; a later `git` or `gh` that is only an
argument): `grep git file`, `grep -rn git .`, `grep -c gh tools/guard.py`,
`rg -n gh tools/`, `echo gh`, `echo git`, `which git gh`,
`brew upgrade git gh`, `git grep -n git -- tools`, `git log --author git`,
`git commit -m "git push is guarded"`, `grep -rn "git push" .claude/`,
`rg 'gh pr merge' tools/`, `grep -rn sh .`, `grep bash -c x`; and reads
through a runner: `caffeinate -i git status`, `caffeinate -i gh pr view 12`,
`watch -n5 'gh pr checks 185'`, `watch -n5 gh pr checks 185`.

Pinned false positives (ask, PUSH): `echo git push`, `grep git push file`,
`man git push`, `parallel 'echo git push' ::: a`. An ask costs a click by
day; these words are rare unquoted.

Pinned passes (residual, §6): `caffeinate -i git p origin` and
`caffeinate -i gh pm 12` (an alias through a bare tail; rule a drops G2's
reason), and `python3 -c "import subprocess; subprocess.run(['git',
'push'])"` (a program, not an argv; it passes on `28b8292` too).

**Red test, review round 8** (`test_guard_push.py`, both modes; each
passes on `fa9f3e1` and asked in this design's in-process prototype of the
two amendments, run against `fa9f3e1`'s `runs`, `publishes` and
`segment_why`, not through the hook; the prototype was not committed):

- [PUSH] `dash -c 'git push'`, `ksh -c 'git push'`, `csh -c 'git push'`,
  `tcsh -c 'git push'`, `fish -c 'git push'`, `/bin/ksh -c 'cd x && git
  push'`, `caffeinate dash -c 'git push'`, `find . -exec tcsh -c 'git
  push' \;` (written as a raw string), `watch "dash -c 'git push'"`;
- [UNKNOWN] `parallel git ::: push`, `parallel ::: git ::: push`,
  `parallel git {} ::: push`, `parallel -j2 git ::: push`;
- a test that every basename listed in the running host's `/etc/shells`
  that ends in `sh` is in `shell_scan.SHELLS`, skipped with a reason where
  `/etc/shells` is absent (read from `/etc/shells`, not `ls /bin/*sh`: on a
  merged-`/usr` Linux runner `/bin/*sh` also lists `ssh`; Ubuntu's bash
  package lists `rbash` in `/etc/shells`, so the set holds `rbash` too;
  not checked here, having no Linux host: CI's run of this test checks
  it); and that
  `guard_push` holds no shell list of its own (`guard_push.SHELLS`, if it
  exists, is `shell_scan.SHELLS`).

Passes (controls): `grep -rn dash .`, `grep ksh file`, `echo tcsh`,
`which dash ksh`, `ls /bin/*sh`, `man csh`, `dash -c 'git status'`,
`parallel 'git status' ::: a`; and `parallel git push ::: a` still asks
with the PUSH reason. Pinned false positives (ask, UNKNOWN): `parallel git
::: status`, `parallel grep git ::: a`.

`test_guard_governance.py`: `dash -c 'cp x CLAUDE.md'` and `ksh -c 'echo x
> CLAUDE.md'` ask (both pass on `fa9f3e1`).

The `\;` strings: two tests on `fa9f3e1`
(`tests/python/test_guard_push.py@fa9f3e1:627` and `:644`) hold `"\;"` in a
plain string, which Python 3.14 warns about. `@tester` makes them raw
strings in this red commit: the green commit touches no test file
(`docs/increments/README.md`, step 3), so review round 8's "in the same
green" is taken by the red step instead.

**Size.** About 15 lines; after review round 8 about +1 (the shell set
grows by one line, `PROGRAM_FLAG` is derived from it, `guard_push.py`'s
copy is dropped, and rule (a)'s `bare` flag changes in place).

**Amendment after review round 9** (Ola's option A, last section).
Rule (a)'s `parallel` exception and rule (c) look only at the command's own
program, so a runner put in front of `parallel`, `watch` or `flock` hides
it: on `a0001d00`, `caffeinate parallel git ::: push`, `caffeinate watch
'git push'`, `stdbuf -o0 flock /tmp/l -c 'git push'` and every other row
of the matrix below with a runner that is not a wrapper pass (36 of its
96 rows, judged with the guard's own `runs`, `publishes` and
`segment_why`; 20 of them also through the hook). Wrappers that `shell_scan` strips
(`nohup`, `time`, `env` and the rest of `WRAPPERS`) are not affected; they
ask today.

The same function has a second defect, found while designing this
amendment: rule (b) parses the rest of the line at every shell word, and
the parsed command is itself scanned for shell words, so a line of *n*
shell words costs about 2^*n* parses. On `a0001d00`, through the hook, 23
`sh` words before `-c 'git push'` took 6.4 s in one run and 50.6 s in
another, and 40 did not finish in 120 s; the time doubles with each word.
Claude Code's `PreToolUse` command hooks time out at 600 s, and "a
timed-out `command` … hook doesn't block the tool call"
(code.claude.com/docs/en/hooks, read 2026-10-06), so about 30 `sh` words
before a push pass the guard unasked. On `master` (`ed125121`) the same
line takes 0.07 s and passes, since rule (b) is new in this PR. Rule (d)
below needs the same shape as the fix, so the fix costs no extra line;
§7 question 5 asks Ola whether it stays in this round (default: yes).

**Change.** `runs` takes two flags that pass down the chain of programs:

    def runs(words: list[str], under: bool = False, quoted: bool = False)
            -> list[tuple[list[str], bool]]
        # under: this command, or a runner before it, is `parallel`;
        # quoted: this command, or a runner before it, is in STRING_RUNNERS.

Its first lines set `under |= name == "parallel"` and `quoted |= name in
STRING_RUNNERS`. Then, for each word after the first:

- **(a)** unchanged but for the flag: a bare `git` or `gh` tail has `bare`
  equal to `not under`.
- **(c)** keyed on `quoted`, not on the command's own program: each later
  word with whitespace is parsed and each command it returns goes through
  `runs` with no flags (a command string starts a new chain).
- **(d), new: a later word whose basename is in `STRING_RUNNERS`.** The
  rest of the line from it is a command of its own:
  `runs(words[at:], under, quoted)`, and the scan stops there, since that
  call reads every later word. So `caffeinate parallel git ::: push` is
  judged as `parallel git ::: push` (unknown command), and `caffeinate
  watch 'git push'` as `watch 'git push'` (push).
- **(b), reshaped the same way:** at a shell word, the rest of the line is
  `runs(words[at:], under, quoted)`, and each command of the shell's own
  program (the commands after the first that `shell_scan.parse` returns
  for `shlex.join(words[at:])`) goes through `runs` with no flags; then
  the scan stops. A word is a shell or a string runner, never both, so (b)
  and (d) do not meet.

The depth of the recursion is the number of runner and shell words on
the line. The time is not linear (corrected after review round 10): each
shell word parses the rest of the line again (`shell_scan.parse` in (b)),
so the time grows with the number of shell words times the length of the
line. Measured by `@reviewer` through the hook at `80800a14`: 400 `sh`
words before `-c 'git push'`, 0.16 s; 800, 0.39 s; 800 followed by 100 KB
of plain words, 31 s; the same line with no `sh` words, 0.08 s. Its
cause, the slow line of §7 question 5, grew with the number of shell
words alone (about 30 of them passed 600 s); this one needs a long line
as well. **Constant:** Python's default recursion limit, 1000
frames. The deny starts at about 990 runner or shell words; the exact
count depends on the form of the line. Measured through the hook at
`85d39f5c` with Python 3.14.7: 990 `sh` words before `-c 'git push'`;
989 `caffeinate sh -c` groups before `'git push'`; 993 `watch` words
before `'git push'` (the command quoted), but 999 before `git push`
unquoted, which at 993 to 998 still asks. Review round 11 measured 995
for the quoted `watch` form, so a count can differ by a few between
setups. Past the limit, the call raises `RecursionError`, which the
guard's `except Exception` turns into a deny (`guard_push failed:
RecursionError: …`), so it fails closed, but with a reason worded as a crash. A cap on
runner and shell words per line would bound the time and give a plain
reason; it is in §6, for a later harness increment.

**Why not memoise instead.** Caching `runs` by its words keeps the parse
count linear but not the result: each call's list holds every later
call's, so the list itself doubles per word. Stopping the scan is what
bounds the list; the parse count is bounded only by the line (above).

Probed with a scratch prototype of the reshaped `runs` (about +3 net
lines on `guard_push.py`), in a scratch copy of `a0001d00` made by
`tools/scratch_copy.py`, removed: the five guard suites passed (1046
tests), all 96 rows of the matrix below asked with the head's reason (60
of them also through the hook), the time rows answered in under 0.2 s
each, and the controls and pins behaved as listed.

**Red test, review round 9** (`test_guard_push.py`, both modes):

- **The runner matrix.** Each runner placed before each head asks with
  the head's reason. Runners: `caffeinate`, `caffeinate -i`, `stdbuf
  -o0`, `arch -arm64`, `find . -maxdepth 0 -exec … \;` (a raw string),
  `env FOO=1 caffeinate`, and every name in `shell_scan.WRAPPERS`, read
  from the module and not copied (`timeout 5`, `script -q /dev/null` and
  `sudo` with their arguments; the wrapper rows ask on `a0001d00` too and
  guard the stripping). Heads: `parallel git ::: push` [UNKNOWN],
  `parallel ::: git ::: push` [UNKNOWN], `parallel gh ::: pr ::: merge`
  [UNKNOWN], `parallel 'git push' ::: a` [PUSH], `watch 'git push'`
  [PUSH], `flock /tmp/l -c 'git push'` [PUSH]. Red: the six non-wrapper
  runners on every head.
- [UNKNOWN] `watch parallel git ::: push`, `caffeinate watch caffeinate
  parallel git ::: push`; [PUSH] `parallel sh -c {} ::: 'git push'`
  (asks on `a0001d00`; keeps (c) alive through a shell under the new
  flags).
- **Time, measured from outside the hook process** (a subprocess with a
  timeout): `sh ` forty times then `-c 'git push'`, `watch ` forty times
  then `'git push'`, and `caffeinate sh -c ` thirty times then `'git
  push'`, each answered `ask` [PUSH] within 5 s; `bash ` four hundred
  times then `-c 'git status'` passes within 5 s (no deny from the
  recursion limit). The first fails on `a0001d00` by timeout.

Passes (controls): `caffeinate parallel 'git status' ::: a`, `stdbuf -o0
flock /tmp/l -c 'git status'`, `caffeinate watch -n5 'gh pr checks 185'`,
`caffeinate watch -n5 git status`, `watch sh -c 'git commit -m "git
push"'`, `git log --grep 'watch git push'`; round 7's and round 8's
controls stay as they are.

Pinned false positives (ask): `grep parallel git ::: x` and `echo
parallel gh ::: pr` [UNKNOWN], `caffeinate parallel git ::: status`
[UNKNOWN, as `parallel git ::: status` already is], and `grep watch 'git
push' file` [PUSH]: a runner's name as an argument now starts a chain.

**Size, round 9.** About +3 (two flag lines, (d) one line, (b) reshaped in
place).

**The governance half stays in §6.** `guard_governance.py` would need the
same tail rule with every writer as a tail start (`cp`, `mv`, `tee`, `ln`,
`install`, `rm`, `touch`, `truncate`, `sed -i`, `dd`), and those are common
arguments: `writer_targets('rm', ['CLAUDE.md'])` returns `['CLAUDE.md']`
(checked on `28b8292`), so `grep -n rm CLAUDE.md` would ask, and
`grep tee .claude/agents/x.md` likewise, each a refusal while unattended.
Rule (a)'s two names do not have that problem. The governance Bash arm is
also already open to any program that computes a path (G5), and G5's brief
line moves rule-file writes to Edit and Write, which judge the exact path;
a runner opens no route of a new kind there. A fix belongs in
`tools/shell_scan.py`'s `unwrap`, with each runner's option grammar, after
the freeze. (The widened shell set does reach this guard, since a shell
word is parsed as a nested line, not unwrapped: rule (b)'s red cases.)

### G8. A copy or move into a directory named without `/` (review round 8)

**Change.** `tools/shell_scan.py`'s `writer_targets`, for `cp`, `ln`,
`install` and `mv` (now handled in the same branch; `mv` also keeps its
sources, which it removes):

- **A target directory option**, `-t DIR`, `-tDIR`, `--target-directory
  DIR` or `--target-directory=DIR`: every operand is a source, and each is
  written as `DIR/<its basename>`. Today the option's value is skipped as
  an option argument, so `cp -t .claude/hooks x.py` names only `x.py`.
  These are GNU options; this Mac's `/bin/cp` and `/bin/mv` are BSD's and
  have none of them, but the parser does not know which `cp` runs.
- **A last operand with no file extension** (`os.path.splitext` gives
  `''`: `tools`, `.claude/hooks`, `.git/remotes`, `.claude`, `.git`): it is
  judged both as itself and as a directory, `<it>/<basename of each
  source>`. Whether it is a directory is not known before the line runs
  (a `cd` or `mkdir` earlier on the line changes it), so the guard does
  not look; it judges both readings. A last operand with an extension
  (`/tmp/x.md`, `notes.bak`) is read as a file, as today.
- `.`, `..` and a trailing `/` are unchanged (the directory reading only).

**Why here, not in the guard.** Every caller of `writer_targets` gets the
same list of targets; `guard_governance.py` judges each with `governed`, so
G4's `tools/<stdlib name>` rule, G1's `.git/remotes/` and `.git/branches/`
prefixes, every other prefix and the settings glob all apply unchanged.

**Incident.** Review round 8, blocking item 4: `cp json.py tools`, `mv
ast.py tools`, `cp s .git/remotes` and `cp x.py .claude/hooks` pass;
`writer_targets` has judged only the target's own name since h4.

**Red test** (`test_guard_governance.py`, through `judge_bash` or the hook
in both modes; each passes on `fa9f3e1`, checked with `fa9f3e1`'s
`writer_targets` and `governed`, and asked, naming the path in brackets, in
this design's in-process prototype). Asks: `cp json.py tools`
[`tools/json.py`], `mv ast.py tools` [`tools/ast.py`], `cp a.py b/json.py
tools` [`tools/json.py`], `cp s .git/remotes`, `cp s .git/branches`, `cp
x.py .claude/hooks`, `cp a.md .claude/agents`, `cp settings.json .claude`
[`.claude/settings.json`], `cp config .git` [`.git/config`], `ln -s
/x/y.py .claude/hooks` [`.claude/hooks/y.py`], `install x.py
.claude/hooks`, `cp -t .claude/hooks x.py`, `cp -t.claude/hooks x.py`,
`cp --target-directory=.claude/hooks x.py`, `mv -t tools ast.py`. A
`shell_scan` unit test: `writer_targets('mv', ['a', 'b', 'tools'])` holds
`a`, `b`, `tools`, `tools/a` and `tools/b`.

Passes (controls): `cp CLAUDE.md /tmp/x.md`, `cp notes.txt notes.bak`,
`cp json.py src`, `cp x.py tools`, `mv ast.py tools/ast_helpers.py`, `mv
tools/old.py tools/new.py`, `cp -r docs/increments /tmp/inc`; and `cp
notes.txt CLAUDE.md` and `cp x.py .claude/hooks/` still ask.

Pinned false positive (asks): `cp CLAUDE.md backup`, where `backup` is
meant as a new file: read as a directory, it would write
`backup/CLAUDE.md`, and a file named `CLAUDE.md` is governed anywhere.

**Size.** About 8 lines.

### G4. Guard import shadowing

**Change.** Two parts, both needed:

1. Every `sys.path.insert(0, <tools>)` becomes `sys.path.append(<tools>)`:
   `.claude/hooks/guard_push.py@bc01cd8:32`,
   `.claude/hooks/guard_governance.py@bc01cd8:44`,
   `.claude/hooks/guard_spawn.py@bc01cd8:80`,
   `.claude/hooks/guard_unattended.py@bc01cd8:76`,
   `tools/session_state.py@bc01cd8:48`, `tools/brief.py@bc01cd8:28`,
   `tools/away.py@bc01cd8:32` (seven lines changed, net 0). The standard library
   then wins over any file in `tools/`. The hooks' own directory,
   `.claude/hooks/`, stays first, and is governed by prefix.
2. A script run as `python3 tools/x.py` (the `SessionStart` hook, the main
   session's `brief.py`) has `tools/` first by Python's own rule, which no
   `append` changes. So `governed()` also returns True for a path whose
   component after a `tools` component, cut at the first `.`, is in
   `sys.stdlib_module_names` (`tools/ast.py`, `tools/json/__init__.py`,
   `tools/subprocess.cpython-314-darwin.so`). G3a's exemption comes first,
   so a probe in a scratch copy is not asked about.

**Incident.** h12 PR A review round 2 (recorded in 6403582 on
`worktree-h12-design`): a planted `tools/json.py` replaced the stdlib module
inside `tools/ci_changes.py`, fixed there with `python3 -I`. The reviewer's
lesson carried it to the gates; `@architect` then found it reaches the
hooks (main session transcript `806b4380`, line 7527, the `ASK OLA:` line
written at 07:00:27). Checked again on this branch's base in a scratch copy
of the two guards and `tools/`: with `tools/dataclasses.py` raising at
import, `guard_governance.py` given a Write of `CLAUDE.md` and
`guard_push.py` given `git push` both exited 1 with a traceback and no
decision, which Claude Code treats as a non-blocking error; with
`sys.path.append` in both, each printed its `ask`. No file in `tools/` or
`.claude/hooks/` is named after a stdlib module on any branch (the main
session's check, same transcript, line 7548).

**Red test.** `test_guard_governance.py`: `governed()` is True for
`tools/ast.py`, `tools/json/__init__.py`, `/abs/wt/tools/typing.py`; False
for `tools/scratch_copy.py` and `docs/x/ast.py` (`tools/count_loc.py`
until review round 9: it is governed by §3, and the suite already uses
`tools/scratch_copy.py`). A test per hook that copies
the hook and `tools/` into `tmp_path`, plants `tools/dataclasses.py` raising
at import, feeds a governed Write (or a push) and requires an `ask` on
stdout. A repository test: no file in `tools/` is named after a stdlib
module, so one approved by mistake still fails CI.

**Size.** About 6 lines net.

**Amendment after review round 9** (Ola's option A, last section). Part
2 reads the path as written, so `cp json.py tools/.`, `cp json.py
tools/./` and `mv json.py tools/./json.py` pass on `a0001d00`
(`writer_targets` names `tools/./json.py`, and the component before
`json.py` is `.`, not `tools`). The same spelling opens the other rules:
`echo x > .claude/./agents/x.md`, `cp a.md .claude/./agents/a.md` and a
Write of `<worktree>/.claude/./agents/x.md` or `<worktree>/tools/./json.py`
pass there too (each probed on `a0001d00`, the Write through the hook).

**Change.** `governed()` judges the path also in its normalised form:
after computing `norm`, if `posixpath.normpath(norm)` differs from `norm`
and is not `.`, `governed(posixpath.normpath(norm), scratch_exempt=…)`
returning True makes the path governed. Both forms, not the normalised
one alone: `normpath` drops a trailing `/` and `/.`, so `rm -r
.claude/hooks/` and `rm -r .claude/hooks/.`, which ask today by the
`.claude/hooks/` prefix, would pass if only the normalised form were
judged. G3a's scratchpad check stays first. `normpath` is lexical: a
`..` after a symbolic link is resolved as text, not as the system would
(§6). About +3 lines (an import and the check).

**Red test, review round 9** (`test_guard_governance.py`; each passes on
`a0001d00` and asked in a scratch-copy prototype of this change). Asks:
`cp json.py tools/.`, `cp json.py tools/./`, `mv json.py
tools/./json.py`, `cp json.py tools//json.py`, `cp json.py
tools/x/../json.py`, `cp json.py ./tools/.`, `cp -t tools/. json.py`,
`install json.py tools/.`, `echo x > tools/./ast.py`, `echo x >
.claude/./agents/x.md`, `cp a.md .claude/./agents/a.md`; a Write of
`<worktree>/.claude/./agents/x.md` and of `<worktree>/tools/./json.py`
through the hook, both modes. `governed()` is True for `tools/./json.py`,
`tools//ast.py`, `a/tools/x/../typing.py`, and False for
`tools/./scratch_copy.py`, `docs/x/./ast.py`, `tools/../json.py`. Still ask:
`rm -r .claude/hooks/` and `rm -r .claude/hooks/.`. Passes (controls):
`cp x.py tools/.`, `cp json.py src/.`, `mv tools/./old.py
tools/./new.py`, `cp json.py tools/../json.py`, `cp CLAUDE.md /tmp/x.md`,
and a Write of `<scratchpad>/tools/./json.py` (G3a).

### G5. Rule files written through the shell

**Finding: the reported bypass did not happen.** The h12 `@architect`
(agent `a329f25e`, the run that committed 6403582) handed back "I wrote two
governed files with a Python heredoc, which the governance guard does not
see." Its transcript shows otherwise: the heredoc that wrote
`docs/increments/README.md` and `.claude/agents/reviewer.md` (line 344,
06:49:26 UTC) drew `guard_governance.py`'s `ask` naming both files (line
346), which was answered yes 71 seconds later (line 348); its reverts were
asked about too (lines 355, 364). Replayed on this branch's base,
`judge_bash` on that command returns `ask` for both files.

**What is true.** The Bash arm reads string literals, so a path built from
parts passes: on this branch's base, `judge_bash` returns None for
`Path('docs', 'increments', 'README.md').write_text(…)`,
`(Path('docs/increments') / 'README.md').write_text(…)`,
`open('CLAU' + 'DE.md', 'w')`, and `python3 /tmp/w.py`. The guard's own
docstring calls itself "a tripwire, not a sandbox"; no static reading of a
program closes this.

**Change (bound, not closed).** One line in `.claude/briefs/common.md`:
"Write a rule file with Edit or Write, never through the shell." The Edit and
Write arm judges the exact path, so the rule moves rule-file writes to the
arm that cannot miss. Prose, no red test. (G5b, joining literal path parts
in `tools/shell_scan.py`'s `candidates`, about 15 lines, is left out, §5.)

### T1. `tools/count_loc.py` (N1)

**Interface.**

    python3 tools/count_loc.py <base> [<head>]

`<head>` defaults to `HEAD`. The old side is `git merge-base <base> <head>`,
so `<base>` may be `origin/master` and the count is the PR's, as
`git diff <base>...<head>`. Output, tab-separated, one line per counted
file, then a total, then one `not counted:` line per changed file it skipped
with the reason (`tests/`, `docs/`, `not code`):

    src_python/tin_engine/gauge.py	117	0	117
    total	744	53	691
    not counted: tests/python/test_gauge.py (tests/)

Exit 0; exit 2 with one stderr line starting `count_loc:` on a git error.
Committed revisions only; the working tree is not read.

**Rule, from `CLAUDE.md` §2.** Added lines are the `+` ranges of
`git diff -U0` hunks, judged in the file at `<head>`; removed lines the `-`
ranges, judged at the old side. A line counts unless it is blank, a comment,
a docstring, or inside a raw literal's body. Per kind:

- Python (`.py`): `tokenize`, as the 29 PR 2 round-2 record states its
  method (cc52b8f): a line counts if it holds a token other than comments,
  newlines and indentation; docstrings (the first statement of a module,
  class or function, if a string; found with `ast`) do not count on any of
  their lines; any other string spanning lines counts on its first line
  only. A file that does not tokenize counts every non-blank line, with a
  stderr warning.
- C++ (`.h`, `.hpp`, `.cpp`, `.cc`, `.cxx`): a small scanner; a line counts
  if it has a character outside `//` and `/* */` comments and outside a raw
  string's body. A raw string `R"d(…)d"` counts on its opening line; its
  later lines count only for code after the closing `)d"`.
- CMake (`CMakeLists.txt`, `.cmake`) and shell (`.sh`): blank and
  `#`-comment lines do not count.
- Not counted: anything under `tests/` or `docs/`, and files of no kind
  above (Markdown, YAML, TOML, JSON, data). The `not counted:` lines show
  them, so a reviewer sees what was left out (question 2).

Every git call runs with `GIT_CONFIG_GLOBAL=/dev/null`,
`GIT_CONFIG_NOSYSTEM=1`, `GIT_NO_REPLACE_OBJECTS=1`, `--no-ext-diff
--no-color`: a user setting such as `diff.noprefix` changes the headers the
parser reads (h12 §3.3 names it). Renames are git's default detection.
The repository's own config is not neutralised by those variables, so the
diff also passes its settings explicitly: `--src-prefix=a/ --dst-prefix=b/`
and `--inter-hunk-context=0` (review round 1), and `--no-relative` (review
round 2: `diff.relative=true` in `.git/config`, run from a subdirectory,
would otherwise count only that subdirectory).

**Blueprint.** Pure functions, each tested alone, and one thin `main`:

    def kind_of(path: str) -> Kind | None
    def counted_lines(text: str, kind: Kind) -> frozenset[int]   # 1-based
    def hunks(diff: str) -> dict[str, tuple[list[int], list[int]]]  # path -> (new +lines, old -lines)
    def tally(...) -> list[Row]                                    # joins the three
    def main(argv) -> int                                          # git I/O only here

**Incident.** Night retrospective, lesson 3 (a79f2d9): eighteen runs wrote
their own counter; two personas got 558 and 577 from the same tree
(`a1fee960`'s lessons).

**Red test** (`tests/python/test_count_loc.py`): `counted_lines` cases per
kind (comment-only, blank, docstring of one and of three lines, a
multi-line non-docstring string, a C++ raw string over three lines with
`;` after it, `/* */` across lines with code after it); `hunks` on a
hand-written diff with a new file, a deletion and two hunks; a temporary
repository whose two commits give known counts; that a `GIT_CONFIG_GLOBAL`
file with `diff.noprefix=true` does not change the output; and two
recorded counts from real history, which need the full clone CI already
has (`fetch-depth: 0`): `count_loc.py 529613a 193079d` totals 589, 20,
569 (29 PR 2, review round 3), and `count_loc.py 9e666f4 9bb1723` totals
744, 53, 691 (29 PR 4, ef403bd's message and the round-6 review). If the
tool and a recorded count disagree, `@tester` reports which line differs
before anyone changes either.

**Size.** About 120 lines.

### T2. `tools/scratch_copy.py` (N2), and P9

**Interface.**

    python3 tools/scratch_copy.py <rev> <dir>

1. Refuses (exit 2, one stderr line) if `<dir>` exists and is not empty, or
   lies inside the worktree it is run from or inside the main checkout (the
   work tree of the common git dir, which holds `.claude/worktrees/`).
2. `git archive <rev> | tar -x -C <dir>`, with T1's three variables, so the
   copy is not a git work tree.
3. Copies the built `_core*.so` from the running worktree's
   `.venv/lib/python3.*/site-packages/tin_engine/` into
   `<dir>/src_python/tin_engine/`. None found: a stderr warning (Python-only
   suites still run). If `git diff --quiet <rev> HEAD -- include src
   bindings CMakeLists.txt` says the C++ differs, a stderr warning that the
   copied `_core` was built from `HEAD`.
4. Writes `<dir>/.scratch_copy/sitecustomize.py`, which drops the editable
   finder from `sys.meta_path`, and prints one line on stdout, the command
   that runs `pytest` against the copy: `cd <dir> && PYTHONPATH=<dir>/.scratch_copy:<dir>/src_python`
   and the worktree's `.venv/bin/python -c` with a program that only calls
   `pytest.main` on the remaining arguments (`tests/python/` by default;
   replace it to run some tests only). Because `PYTHONPATH` is inherited,
   the drop runs at every interpreter start, after the `.pth` file installed
   the finder, so `tin_engine` comes from the copy in pytest's process and
   in any Python child a test starts. Review round 1 caused this: the first
   version did the drop inside the `-c` program only, as
   `docs/benchmarks/2026-10-02/basin-memory-probe/blockprobe.py:23` does,
   and a child process reloaded the finder and imported the worktree's
   code. Limits: a child whose environment replaces `PYTHONPATH` (one built
   from scratch, as `tests/python/test_io_geotiff.py@c3817b3:1271` does) or that
   starts Python with `-I` or `-E` (both ignore `PYTHONPATH`) still imports
   the worktree's code, so a mutant run there can report a false survivor.

It also serves "run the new tests against the code before the change"
without `git stash`.

**P9** (evening retrospective, c0e3451): `test_hook_is_executable_in_the_checkout`
in `tests/python/test_settings_wiring.py` skips, with the reason stated,
when the tree is not a git work tree. `@tester`'s, in the red step.

**Incident.** Night retrospective lessons 1 and 2 and §3 (a79f2d9): six
`@tester` runs each built their own way round the editable finder; the
seven `test_settings_wiring` failures in a `git archive` copy (round 2's
"8 failed", c7d427a); five `git stash` uses in the shared stash.

**Red test** (`tests/python/test_scratch_copy.py`): the refusals; the copy
has no `.git`; the printed command, run with `tests/python/` replaced by
one generated test file asserting `tin_engine.__file__` lies under
`<dir>`, exits 0 (this is the check that the finder is dropped; it fails
if the drop is removed); the C++-differs warning on a revision before a
C++ change.

**Size.** About 60 lines.

### T3. `tools/brief.py`

**Note-file names.** Today `<persona>-<HHMMSS>.md`. Two briefs printed in
the same second for the same persona in different worktrees got the same
name: `architect-094230.md`, for `h15-ci` and for `h16`, on 2026-10-05
(main session transcript `912df417`, lines 147 and 155; no commit). The
name becomes `<persona>-<worktree>-<HHMMSS>.md`, `<worktree>` the
worktree directory's name; if that file exists, `-2`, `-3`, … are appended
before `.md`. `.claude/REQUIRED-READING.md` states the format and changes
with it. **Red test** (`test_brief.py`): two briefs with the clock frozen,
for `h15-ci` and `h16`, name different files; with the file already there,
the next name ends `-2.md`. About 6 lines.

**Root Markdown files.** The item is the evening retrospective's P1 (c0e3451): root Markdown files
(`README.md`, `INSTALL.md`, `NOTICE.md`, `testing.md`,
`project_structure.md`, `auto_catchments.md`, `parallel_refinement.md`)
are in no persona's write limit, and Ola ruled on 2026-10-05 that
`@architect` may edit them (recorded in ef403bd's message). `WRITES` still
lacks them: this run's own brief prints `@architect`'s limit without them.
The change: `@architect`'s entry gains "root *.md files". P1's second half,
a test that every tracked path falls in some limit, needs the limits as
patterns rather than prose, about 30 more lines; it is left out (§5).
Question 1 asks Ola whether this is the item. **Red test**: `brief.py
--persona architect` prints a limit naming root Markdown files. About 1 line.

### R1. Rule lines

- **Monitor, not `sleep`.** `.claude/briefs/common.md` gains: "Wait for a
  background run with the Monitor tool, never with `sleep`." Incident: night
  retrospective §2 (a79f2d9), three refusals by Claude Code's own `sleep`
  check (`ac6626cc` 00:48:08, `a9fbff4d` 01:54:49, `a79b3df8` 02:50:22),
  each an agent sleeping and then reading its log.
- **A fallback for the last hour.** `.claude/REQUIRED-READING.md`,
  *Unattended mode*, already requires a no-ruling fallback in `QUEUE:`
  before Ola leaves; it gains "and one such item stays queued for the
  window's last hour". Incident: night retrospective §1 and question 1
  (a79f2d9): from about 04:35 to 05:28 every queued item waited on Ola.
- **Pointers.** `tester.md`'s mutant paragraph replaces its `git archive`
  recipe and place list with "in a scratch copy made by
  `python3 tools/scratch_copy.py <rev> <dir>` in the session scratchpad,
  removed afterwards" (the night's cut C4, about 35 words fewer);
  `reviewer.md`'s LOC item and `CLAUDE.md` §2 name `tools/count_loc.py`;
  REQUIRED-READING's *The harness* names G1 to G4 in one sentence (PR B),
  rewritten after Ola's option C to name G1's amendment and G6 and to say
  that G3a holds only for Edit, Write and one plain shell command, and that
  no git write in a scratchpad repository passes; and again after review
  round 7, to name the include and transport overrides and the glued
  `gh api`/`curl` forms, and to say that a git or gh command run by
  another program asks too (G7); and again after review round 8, to give
  the shell example as any shell (`dash -c '…'`) and to say that a copy or
  move into a governed directory asks with or without the trailing `/`
  (G8); and again after review round 9, to say that a runner in front of
  `parallel`, `watch` or `flock` does not hide it (G7), that a path is
  judged also with `.`, `..` and a doubled `/` resolved (G4), and that a
  write to a whole governed directory is not asked about (§6).

Prose; no red test. `@architect` writes them, in the PR that ships the
tool or guard they describe. `CLAUDE.md` changes, so the main session
restarts after PR A merges.

## 3. Shared pieces and boundaries

- **`tools/scratchpad.py`**, new, standard library only: the pattern and
  `def under(path: str) -> bool` (absolute paths only, by `realpath`).
  Imported by `guard_governance.py` (G3a; `guard_push.py` imported it for
  G3b until option C dropped it) and in `GOVERNED`'s self-protecting set,
  as `shell_scan.py` is. About 10 lines.
- **`tools/count_loc.py`** joins `GOVERNED` too: it computes the arithmetic
  of a rule in `CLAUDE.md` §2, so a change to it is a change to the rule.
  `tools/scratch_copy.py` does not.
- The hooks stay pure apart from the git calls named above; the new git
  calls are G2's two command listings, with a fixed argv (G3b's
  `rev-parse` lookup goes with G3b). The guards' `except Exception` → `deny` stays, so a
  failing lookup refuses rather than passes.

## 4. PR split and size

Under `CLAUDE.md` §2 (700 net per PR) everything fits one PR (about 250
lines, about 300 with G3c and G5b). It is split anyway, so the tools are
not held up by guard review rounds (h12's design took four):

| PR | Items | Production lines, about |
|---|---|---|
| A, tools | T1, T2 (+P9 test), T3, G5's line, R1 but its *The harness* sentence | 187 (count_loc 120, scratch_copy 60, brief.py 7) |
| B, guards | G1, G2, G3a, G4, G6, `scratchpad.py`, the `GOVERNED` entries, R1's *The harness* sentence | 60 as first estimated (G1 12, G2 14, G3a 6, G3b 20, G4 6, scratchpad 10, minus shared lines); measured 78 at `084a6b3` (`tools/count_loc.py origin/master 084a6b3`, `origin/master` at `44fa7f5`, PR A merged); after option C about 70 (78, G3b's removal −26, G1's amendment +8, G3a's plain-command test +3, G6 +8); measured 77 at `a42d864` (`tools/count_loc.py origin/master a42d864`, merge base `7dde17a`); after review round 7 about 88 (77, G1's include/transport keys, text check and governed prefixes +3, G6's glued `gh api`/`curl` options +8); with G7 (Ola's ruling on §7 question 4) about 103 (+15); measured 97 at `fa9f3e1` (`tools/count_loc.py origin/master fa9f3e1`, merge base `7dde17a`); after review round 8 about 106 (G7's shells and `parallel` +1, G8 +8); measured 109 at `a0001d00` (`tools/count_loc.py origin/master a0001d00`, merge base `879ea493`); after review round 9 about 115 (G7's runner rule (d) and the reshaped (b) +3, G4's normalised path +3) |

PR B's branch took `origin/master` in before round 8's red step (merge
`f86d3953`), so the guard test files carry the `harness` pytest marker as
CI runs them. `origin/master` has moved on since (`git log --oneline
HEAD..origin/master` lists what the branch lacks); round 9's red and green
steps do not need it, and whether to merge again before the push is the
main session's call.

Order: A first. Each PR runs red (`@tester`), green (`@developer`), review
(`@reviewer`). Neither touches refine or mesh code, so no `@perf` run. No
suite here is invariant-critical, so no mutation round. Writing the hooks,
`tools/` files in `GOVERNED` and the rule files is asked about at each
write, so both PRs are day work.

## 5. What should drop

- **G3c (the push to a scratch repository): drop, by default.** It is the
  one part that lets a push through the guard, its safety rests on the
  route list being complete, and h12's review found a new route in each of
  three rounds. The night's push refusals had a workaround that cost
  seconds (`git remote get-url --push`). Question 3.
- **G5b (joining literal path parts): drop.** The incident did not happen
  as reported, and what remains open after it stays open after it.
- **P1's coverage test: drop** from h16 (above).
- **G3b (local git writes in a scratch repository): dropped** after review
  round 6, by Ola's option C (last section).

Nothing else should drop: G4 is a live hole in the hooks, G1 and G2 are
small, and T1 and T2 remove the most repeated work of the night.

## 6. Residual, after h16

The guards remain tripwires. Not covered: programs that compute a path
(G5); a module planted in the user-writable Homebrew site-packages under a
`tools/` module's name, which with `append` would now win (h12 §5's
environment residual); shell startup files and `GIT_*` variables they
export, including `GIT_CONFIG_*` set by an earlier command (G1's amendment
reads only the line's text); in a scratchpad, a single command that makes a
link and writes through it, and a background job that re-points a link
between the hook's check and the write (G3a); relative paths in a
scratchpad, which are judged as anywhere else.

Added after review round 7, each checked on `a42d864`:

- **Command runners `tools/shell_scan.py` does not unwrap.** It strips
  `env`, `nohup`, `time`, `timeout`, `nice`, `command`, `exec`, `xargs`,
  `script` and `sudo`; any other program that runs its arguments as a
  command hides that command from both guards. Through `guard_push.py`'s
  hook, each of these printed nothing (a pass), though each runs a plain
  `git push`: `find . -maxdepth 0 -exec git push origin HEAD \;`,
  `caffeinate -i git push`, `stdbuf -o0 git push`, `watch -n1 git push`,
  `flock /tmp/l git push`, `parallel git push ::: a`, `arch -arm64 git
  push`. (`git rebase --exec "git push"` and `xargs git push` ask.) The
  same runners hide a governed write from `guard_governance.py`:
  `caffeinate -i cp notes.txt CLAUDE.md`, `stdbuf -o0 cp notes.txt
  CLAUDE.md` and `find . -maxdepth 0 -exec cp notes.txt CLAUDE.md \;` pass,
  where `cp notes.txt CLAUDE.md` asks. Question 4. [Ola ruled to close
  it now: the push half is closed by G7. What stays here: the governance
  half, for the reasons at the end of G7; a git or gh alias run by another
  program (`caffeinate -i git p origin`, `caffeinate -i gh pm 12`), since
  a bare tail does not apply G2's unknown-command reason; a command string
  given to a runner other than `watch`, `parallel` and `flock` that is not
  a shell (`tmux new 'git push'`, `ssh host 'git push'`), unless a later
  bare `git` or `gh` word catches it; and an interpreter program that runs
  git (`python3 -c "import subprocess; subprocess.run(['git', 'push'])"`),
  which passes the push guard as it passes the governance guard (G5).]
- **Config that runs a program**, on any git command, not only a fetch:
  `-c core.fsmonitor=<program>`, `-c core.hooksPath=<dir>` (a
  `reference-transaction` or other hook in it), `core.askPass`,
  `core.alternateRefsCommand`, `uploadpack.packObjectsHook`,
  `credential.helper`, and `core.sshCommand` or `GIT_SSH*` on a command
  other than fetch (fetch now asks on them, G1). Whatever such a program
  does is invisible to the guards.
- **User-level config, not governed:** `~/.gitconfig` and
  `~/.config/git/config` (a `remote.`, `url.` or `include.` entry there
  re-points a plain `git fetch origin`, as `HOME=` does for one command),
  and `~/.curlrc`; and `curl -K <file>`, whose file can hold both the forge
  URL and a body, so the forge host is not on the line and the `curl` check
  never runs.

Added after review round 8, each checked on `fa9f3e1` with the push
guard's own functions (`runs`, `publishes`, `segment_why` over
`shell_scan.parse`), not through the hook; each passes there and on
`a42d864` (the round's finding):

- **The forge check's spelling:** `/usr/bin/curl -XPUT
  https://api.github.com/…` (the check matches the word `curl` exactly)
  and `curl -sd x https://API.GITHUB.COM/x` (the host match is
  case-sensitive).
- **A git word the parser cannot see:** `echo push | xargs git` (after
  unwrapping, the command is a bare `git`; its words come from stdin),
  `G=git; $G push`, and the dashed program `$(git --exec-path)/git-push
  origin`.
- **git running a command string:** `git submodule foreach 'git push'`.
- **A fetch into a tag:** `git fetch origin tag v9` writes `refs/tags/v9`
  (it refuses to overwrite an existing tag).

And from this design's amendments: a shell outside `shell_scan.SHELLS`
(none on Ola's Mac; `pwsh`, `nu` and the like), and an interpreter that
runs a command (`tclsh`, `expect`), which pass as `python3 -c` does (G5);
a directory named with an extension (`cp x.py some.d`), read as a file
(no governed directory is so named: G8).

Added after review round 9, each checked on `a0001d00` with
`guard_governance.judge_bash` (passes there, and the first four on
`master` `ed125121` too):

- **A write to a whole governed directory** (Ola's option A: a later
  harness increment). `rm -rf .claude/hooks`, `rm -r .claude/agents`, `mv
  .claude/agents /tmp/` and `cp -R hooks .claude/` pass: the target is the
  directory's own name, which no governed prefix (`.claude/hooks/`, with
  its slash) or file name matches. A fix judges a target that is a
  governed directory, or a parent of one, as governed; `rm -r
  .claude/hooks/` already asks, by the prefix.
- **G8's target-directory option in other spellings** (GNU tools only;
  this Mac's `/bin/cp`, `/bin/mv` and `/bin/ln` are BSD's and have
  neither): `-t` clustered with other flags (`cp -rt .claude/hooks x.py`,
  `cp -vt .claude/hooks x.py`, `ln -st .claude/hooks x.py`, `cp -rt tools
  json.py`), and `--target=`, an abbreviation GNU accepts for
  `--target-directory=` (`cp --target=.claude/hooks x.py`, `mv
  --target=tools ast.py`).
- **G4's normalised path is lexical**: `posixpath.normpath` resolves
  `x/../` as text, so where `x` is a symbolic link the system writes
  elsewhere than the guard judged.

And the new false positives of G7's round 9 rule (d), pinned in its red
test: a runner's name as a plain argument starts a chain, so `grep
parallel git ::: x`, `echo parallel gh ::: pr` and `grep watch 'git push'
file` ask; `caffeinate parallel git ::: status` asks as `parallel git :::
status` does.

Added after review round 10 (Ola's ruling on its suggestion: a later
harness increment, not this PR):

- **No cap on runner and shell words per line.** G7's `runs` parses the
  rest of the line again at each shell word, so its time grows with the
  shell words times the line's length (31 s for 800 `sh` words and 100 KB
  of plain words, measured by `@reviewer` at `80800a14`), and from about
  990 `sh` words on it ends in a `RecursionError` that the guard turns
  into a deny worded as a crash (§2 G7). A fix asks, with a plain reason,
  on any line with more runner and shell words than a small cap; that
  bounds the time and replaces the crash-worded deny.

## 7. Questions for Ola

1. **Which `brief.py` item did you approve as "a limit on root-file size it
   enforces"?** No source found states a size limit. Default: the
   root-Markdown item in T3 (your 2026-10-05 ruling that `@architect` edits
   root Markdown files, written into `brief.py`'s limits).
2. **Which files does the line counter count?** Default: code files (Python,
   C++, CMake, shell) outside `tests/` and `docs/`; workflow YAML,
   TOML and Markdown are listed as "not counted".
3. **Drop the push half of the scratchpad item (G3c) from h16?** Default:
   yes, drop it; scratch config and remote writes with `git -C` still pass.
   [Superseded by option C (Ola's ruling on review round 6): G3b is
   dropped too, so scratch config and remote writes with `git -C` ask like
   any other.]
4. **Command runners (§6, first item added after review round 7): close
   them in PR B, or leave them in §6?** Each lets a plain `git push` through
   the push guard unasked. A small fix: judge, besides the whole argv, the
   tail of it from any later word whose basename is `git` or `gh` (about 3
   lines in `guard_push.py`; a false positive is an ask on a line such as
   `echo git push`). The governance guard's half needs per-runner handling
   and stays in §6 either way. Default: leave all of it in §6, as review
   round 7 placed it, and take it up after the freeze.
   [Answered: close it now (last section). The push half is G7; the
   governance half stays in §6.]
5. **The slow line (§2 G7, amendment after review round 9): fix it in this
   round?** A line with about 30 shell words before `-c 'git push'` keeps
   the push guard busy past Claude Code's 600-second hook limit, and a
   hook that runs out of time does not stop the command, so the push goes
   through unasked. This PR's own new code causes it; `master` has no such
   delay. The fix is the same change that closes the round 9 runner route,
   at no extra lines. Default: yes, fix it in this round.
   [Answered: yes, the default (*Ola's ruling on §7 question 5*, below).
   Fixed by round 9's green
   step, `80800a14`; the time that remains, shell words times line length,
   is in §2 G7 and §6.]

## Ola's rulings

2026-10-05, Ola, verbatim: "yes, go with all three defaults". So: 1 the brief.py item is `@architect`'s write limit gaining root Markdown files (the earlier phrase "root-file limit" was the main session's shorthand; Ola: "I have no idea what a root-file limit would mean."); 2 the line counter counts Python, C++, CMake and shell outside `tests/` and `docs/`, and lists the rest as not counted; 3 G3c, the push to a scratchpad repository, drops out of h16.

## Review

### Round 1: `@reviewer`, code, early, PR A tools part, `bc01cd8..4db1eab`

`@reviewer`'s verdict, word for word:

**Code review, round 1 (early, tools part of PR A), 2026-10-05.** Range `bc01cd8..4db1eab` (design e9e3c46, rulings 093725b, red 11cee8e, partial green 4db1eab). Verdict: CHANGES REQUESTED.

LOC: 252 net production lines by `tools/count_loc.py bc01cd8 4db1eab`: `tools/count_loc.py` 195 and `tools/scratch_copy.py` 57. The workflow, the tests and the docs are listed as not counted. The 57 was recounted by hand, line by line. The estimate was 187 (count_loc 120, scratch_copy 60, brief.py 7). count_loc is 75 over, because the C++ scanner and the git plumbing were not priced; scratch_copy is 3 under. With brief.py's about 7 lines still to come, the total is about 259 against 700. The design names no split seam, and none is needed.

Packing: 11 `# fmt: skip` regions, 6 in count_loc and 5 in scratch_copy. Together they save 43 lines; ruff format would give 222 and 73. Each region keeps one call's arguments or one literal's entries on two lines instead of one per line, and each stays readable. They are packed for density only; none is needed for the ceiling.

The counter reproduces both recorded counts exactly:
- `529613a 193079d` gives 589, 20, 569.
- `9e666f4 9bb1723` gives 744, 53, 691.

On a range I chose, PR #163 (`45acf22^1..45acf22`), it gives 41, 20, 21. Both C++ files were checked by hand:
- `refine.hpp`: 10 added, 6 removed, 4 net.
- `lattice_mesh.hpp`: 31 added, 14 removed, 17 net.

Test strength: I planted 16 faults, one at a time, in a scratch copy. 14 were caught. Two survived:
- removing the explicit `--src-prefix`/`--dst-prefix`;
- dropping C++ character-literal handling.

`ruff` is clean and the prohibited-dependency gate passes. No red-step scaffolding remains. The workflow change gives the full clone (`fetch-depth: 0`) to the `python` job, which is the only job that runs the whole pytest suite; the pack is 40 MiB.

Citations: `check_citations.py` lists three at-risk citations, and all three are quoted review records of earlier revisions, so they stay.
- `h9-spawn-briefs.md:763` cites `test_brief.py:335` as the file stood at 4ee0328, where that line is a concurrency refusal. That is still true.
- `h11-ci-path-filter.md:451`'s `.github/workflows/main.yaml:306-307` [path written out in full by the main session so the citation gate resolves it] was already off on bc01cd8; the `CI result` name was at line 311 there and is at 314 now.
- `h10-merge-queue.md:183` is unaffected.

Blocking:
1. **[Ola]** There is no CI yet; the branch must be pushed and CI green. The 18 `test_brief.py` failures wait on the refused `tools/brief.py` change.
2. **[now]** `scratch_copy.py`: a child Python process imports the worktree's `tin_engine`, not the copy's. The finder drop lives only in the `-c` program, and a child reloads the editable finder at startup. Probe, with the main checkout's venv in a copy of 4db1eab:
   - parent: `…/scratchpad/copy1/src_python/tin_engine/__init__.py`
   - child: `/Users/skavhaug/projects/rasputin/src_python/tin_engine/__init__.py`

   Three suites spawn such children: `tests/python/test_cli_mesh_geographic.py@97eea35:886`, `tests/python/test_features.py@4db1eab:521` and `tests/python/test_io_geotiff.py@4db1eab:1271`. In a mutant run they would test the original code and report a false survivor. Fix: carry the drop into children. One route is a generated `sitecustomize.py` on `PYTHONPATH` in the printed command; a sitecustomize runs after the `.pth` file has installed the finder. `tests/python/test_io_geotiff.py@4db1eab:1271` replaces the whole environment, so where the drop cannot reach, state the limit in the docstring. Add a test in which a child must see the copy.
3. **[now]** `scratch_copy.py:10-11` says "Append test paths to it, or keep the `tests/python/` it ends with". Appending keeps `tests/python/`, so the whole suite runs as well. It should say to replace `tests/python/`.
4. **[now]** `test_recorded_counts_from_real_history` fails in a `git archive` copy ("not a git repository"). That is the same class of failure P9 fixes for `test_settings_wiring.py`, here in a tool built for such copies. Skip it with a reason when the tree is not a git work tree, using P9's `is_work_tree_top`; keep the failure in a shallow clone.
5. **[now]** The status line of `docs/increments/h16-harness-fixes.md` still says "design … next the red step for PR A". Red is 11cee8e and partial green is 4db1eab. Record the measured 252 lines and the reason for the overrun.

Suggestions:
- **[now]** Strip `GIT_CONFIG_COUNT`/`GIT_CONFIG_PARAMETERS` from the environment the counter passes to git, and pass `--inter-hunk-context=0`. Without that, setting `diff.interHunkContext=50` through `GIT_CONFIG_*` in the environment turns PR #163's 41/20 into 63/42; the net stays 21.
- **[now]** Add a test with repo-level `diff.srcPrefix`/`diff.dstPrefix` set, and a C++ case of `'"'` followed by `/*` inside a string. Those are the two planted faults that survived.
- **[now]** Count `.pyi` as Python (`src_python/tin_engine/_core.pyi` is not counted).
- **[now]** `scratch_copy.py` refuses a target only inside the worktree it runs from. A target inside the main checkout, which contains `.claude/worktrees/`, passes; also refuse there.
- **[now]** A failure in `tar` prints more than one stderr line and leaves a partly filled directory.
- **[now]** `count_loc.py` and `scratch_copy.py` are tracked at 100755, while every other `tools/*.py` file is 100644. No rule decides; pick one.
- **[Ola]** A rename from an uncounted path into a counted one counts 0. Probe: a 40-line `docs/proto.py` renamed to `tools/proto.py` counts 0. `CLAUDE.md` §2's literal text (git's default rename detection) gives the same, so a change is a rule change.
- **[Ola]** `CLAUDE.md` §2 still says only "tests excluded", but the counter also leaves out `docs/` and non-code files (ruling 2). The pending §2 pointer line should state both.

Not pushed; no CI.

Taken in the recording commit: blocking item 5 (the status line).

### Round 2: `@reviewer`, code, early, PR A tools part, `4db1eab..731e2a0`

`@reviewer`'s verdict, word for word:

**Code review, round 2 (early, tools part of PR A), 2026-10-05.** Range `4db1eab..731e2a0` (round 1 recorded 1135e1d, red 3de6e1b, green 64e51a2, citation fix 731e2a0). Verdict: CHANGES REQUESTED.

LOC: 280 net production lines by `tools/count_loc.py bc01cd8 731e2a0`: `tools/count_loc.py` 199 and `tools/scratch_copy.py` 81. This round adds 28 net (43 added, 15 removed): count_loc 8/4/4 and scratch_copy 35/11/24. I recounted the scratch_copy part by hand from `git diff -U0` and got the same figures. The estimate was 187. count_loc is 79 over, from the C++ scanner and the git plumbing. scratch_copy is 21 over, from the sitecustomize, the main-checkout refusal and the cleanup after a failed tar; all three came from review and were not priced. With brief.py's roughly 7 lines still to come, the total is about 287 against 700. The design names no split seam, and none is needed.

Packing: there are now 15 `# fmt: skip` regions, 7 in count_loc and 8 in scratch_copy. Four are new this round:
- the `DIFF` tuple;
- the `git rev-parse` call in `_main_checkout`;
- the `tar` call in `_extract`;
- the `quoted` tuple.

Each one keeps one call's arguments, or one tuple, on two or three lines, and each stays readable. Together the regions save about 55 lines; ruff format would give about 229 and 106. They are packed for density only; none is needed for the ceiling.

Round 1's [now] items, each run:
- Blocking 2 (a child process imports the worktree's code). I repeated round 1's probe using the main checkout's venv and a copy of 731e2a0 made by the script. The pytest process and a child `sys.executable -c` both import `…/rv2/copy1/src_python/tin_engine/__init__.py`, and the child's `_core` is the copy's. In the control, the same program run without the printed `PYTHONPATH` imports `/Users/skavhaug/projects/rasputin/src_python/tin_engine/__init__.py` in both processes. Closed. One limit remains: a child started with `-I` still imports the worktree's code. The only such child in the suite is `tests/python/test_hardening.py@a61e848:247`, which loads `_core` by path and not as the package, so no test is affected.
- Blocking 3: the docstring now says to replace `tests/python/`. Closed.
- Blocking 4: in the copy, the 5 history tests skip with the reason ("… is not the top of a git work tree"). In a `--depth 1` clone of the branch they fail: 2 recorded-count cases and 3 PR #163 cases. Closed.
- Blocking 5: the status line was taken in 1135e1d. It is out of date again (see below).
- Config in the environment and hunk context: closed. Planting the removal of `--inter-hunk-context=0` is caught.
- Prefixes and character literal (round 1's two surviving faults): I planted both again in a scratch clone. The prefix fault fails 1 test and the character-literal fault fails 1 test. Closed.
- `.pyi`: closed; removing it is caught.
- Main checkout: closed; removing `_main_checkout()` from the refusal is caught.
- tar: closed. Removing the cleanup is caught, and so is a stderr longer than one line.
- Modes: closed. Both files are 100644, like every other file under `tools/`.

Test strength this round: I planted 10 faults, one at a time. 9 were caught. One survived: removing the stripping of config given in the environment (`ENV_CONFIG`). All 65 count_loc tests pass, because `--inter-hunk-context=0` alone overrides the key the tests use. The stripping does matter. Run from `src_python/` with `GIT_CONFIG_COUNT=1 GIT_CONFIG_KEY_0=diff.relative GIT_CONFIG_VALUE_0=true`, the counter gives 280 as written and `total 0 0 0` with the stripping removed.

Gates: 111 tests pass in `test_count_loc.py`, `test_scratch_copy.py` and `test_settings_wiring.py`. 347 tests pass and 15 skip in the copy, in the suites that run with the main venv's older `_core`. `ruff check` and `ruff format --check` are clean, and the prohibited-dependency gate passes. No red-step scaffolding remains.

Citations: `check_citations.py` exits 0 and lists 7 at-risk citations, re-read as quotations:
- The 3 that are new on this range all sit in this file's round-1 record, as the code stood at 4db1eab (`test_brief.py:335`, `.github/workflows/main.yaml:306` [path written out in full by `@architect` so the citation gate resolves it], `scratch_copy.py:10`). They stay.
- `h3-unattended-u1.md:759` names its revision (`as of 6a19357`). It stays.
- The other 3 are as in round 1.

The 731e2a0 edit to the round-1 record is accepted. It is needed (the bare path makes the gate exit 1 with "no such file"), it does not change the meaning, and it is bracketed and attributed.

Blocking:
1. **[Ola]** There is no CI yet; the branch must be pushed and CI green. The 18 `test_brief.py` failures wait on the refused `tools/brief.py` change.
2. **[now]** Design §2 T2 no longer describes the code, and `scratch_copy.py`'s docstring points to it as its spec. Step 1 says the script refuses a target only "inside the repository it is run from", but it now also refuses a target inside the main checkout. Step 4 says the command uses a `-c` program that "drops the editable finder … puts `<dir>/src_python` first on `sys.path`". The command now sets `PYTHONPATH` to a generated `sitecustomize.py` and to `<dir>/src_python`, and its `-c` program only calls `pytest.main`. Rewrite steps 1 and 4 to match, and say that review round 1 caused the change.
3. **[now]** The status line still says "partial green `4db1eab` (252 net …)" and "code review round 1 (early) asked for fixes, in progress". It should state the round-1 green (64e51a2), the 280 lines, and review round 2. The recording commit can take this.

Suggestions:
- **[now]** Add a test that fails when the environment-config stripping is removed. One case is `diff.relative=true` given through `GIT_CONFIG_COUNT`, with the counter run from a subdirectory.
- **[now]** The same key in the repository's own config is not neutralised. `diff.relative=true` in `.git/config`, with the counter run from a subdirectory, should count 0 by the same mechanism as above. I could not plant it: unattended mode refused the `git config` write in my scratch clone. Pass `--no-relative` in `DIFF`, or run git from the top level, and add the test.
- **[later]** The generated `sitecustomize.py` hides Homebrew's own `sitecustomize.py`, so in the copy's processes `sys.base_prefix` is not rewritten and the tk path is not added. No test I ran depends on either. Running the hidden file first, or stating the effect in the docstring, would make it explicit.
- **[later]** "In any Python process a test starts" leaves out children started with `-I` or `-E`, which ignore `PYTHONPATH`. The docstring names only children that replace the environment.
- **[later]** When tar fails on a missing target, `_extract` removes the target but keeps any parent folders that `mkdir(parents=True)` created.

Not pushed; no CI.

Taken in the recording commit: blocking item 3 (the status line).

### Round 3: `@reviewer`, code, early, PR A tools part, `731e2a0..6d35fbd`

`@reviewer`'s verdict, word for word:

**Code review, round 3 (early, tools part of PR A), 2026-10-05.** Range `731e2a0..6d35fbd` (round 2 recorded a61e848, design c3817b3, red 9a1f07a, green 6d35fbd). Verdict: CHANGES REQUESTED.

LOC: 280 net production lines by `tools/count_loc.py bc01cd8 6d35fbd`: `tools/count_loc.py` 199 and `tools/scratch_copy.py` 81. This round adds 0 net (1 added, 1 removed): the `DIFF` tuple's second line gains `"--no-relative"`. I checked this by hand against `git diff -U0`; the other 6 changed lines are comments. The estimate was 187. The overrun is as round 2 recorded it. The design names no split seam, and none is needed. No new `# fmt: skip` region.

Round 2's items, each run:
- Blocking 2 (design §2 T2 out of date): closed. Step 1 now names the main checkout. Step 4 names the generated `sitecustomize.py`, the `PYTHONPATH` and the `-c` program that only calls `pytest.main`, and it says review round 1 caused the change. I ran `scratch_copy.py 6d35fbd` into a scratch folder. The printed command and the `.scratch_copy/sitecustomize.py` it writes match step 4 exactly. A target in the main checkout and a target in the worktree are both refused, as step 1 says. With no built `_core` in this worktree, the step-3 warning is printed. The two citations the rewrite added (`tests/python/test_io_geotiff.py@6d35fbd:1271`, `blockprobe.py:23`) point at the subprocess call and at the line that drops the editable finder.
- Blocking 3 (status line): closed in a61e848.
- Suggestion: config in the repository (`diff.relative`): closed. I removed `--no-relative` by hand to check the test catches it, and the `repository config` case fails.
- Suggestion: a test that fails when the environment-config stripping is removed: **not closed**. I removed the stripping again (`_env` keeps every `GIT_CONFIG_*` key) and all 68 count_loc tests pass. The new `GIT_CONFIG_COUNT` and `GIT_CONFIG_PARAMETERS` cases fail only when `--no-relative` is removed as well: with both removed, 3 tests fail. `--no-relative` now overrides `diff.relative` however it is given, so the stripping is again caught by no test. The stripping still matters. With it removed, `GIT_CONFIG_COUNT=1 GIT_CONFIG_KEY_0=diff.algorithm GIT_CONFIG_VALUE_0=histogram` (or `patience`) turns `count_loc.py 9e666f4 9bb1723` from 744/53/691 into 746/55/691. As written, it stays 744/53/691.

Gates: 114 tests pass in `test_count_loc.py`, `test_scratch_copy.py` and `test_settings_wiring.py`. `ruff check` and `ruff format --check` are clean, and the prohibited-dependency gate passes. No red-step scaffolding remains.

Citations: `check_citations.py` exits 0 and lists 10 at-risk citations. 7 are as in round 2. The 3 new ones are on line 631: round 2's record quoting round 1's citations, as the code stood at 4db1eab. They stay.

Blocking:
1. **[Ola]** There is no CI yet; the branch must be pushed and CI green. The 18 `test_brief.py` failures wait on the refused `tools/brief.py` change.
2. **[now]** The comment above `RELATIVE` in `tests/python/test_count_loc.py` says "the environment forms are what the counter's stripping of `GIT_CONFIG_*` neutralises (no other test fails without it)". Since 6d35fbd that is false: `--no-relative` neutralises them too, and removing the stripping fails no test. Correct the comment. Then close round 2's suggestion with a test that fails when the stripping is removed. The direct way is a unit test that `_env()` drops `GIT_CONFIG_COUNT`, `GIT_CONFIG_KEY_n`, `GIT_CONFIG_VALUE_n` and `GIT_CONFIG_PARAMETERS`. Unlike an end-to-end key, it does not stop catching the removal when the next flag is pinned. A `diff.algorithm` case run through the counter is the end-to-end alternative.
3. **[now]** The status line still says "code review round 2 asked for fixes, in progress" and does not name red 9a1f07a and green 6d35fbd. The recording commit can take this.

Suggestions:
- **[now]** The same `diff.algorithm` key in the repository's own `.git/config` is probably not neutralised, because the stripping covers only the environment. In my probe above, added and removed changed and net did not. Passing `--diff-algorithm=myers` in `DIFF`, with a repository-config test, would pin it; myers is git's default, so the recorded counts stay. I could not test this: unattended mode refused the `git config` write in my scratch clone.
- **[later]** Round 2's three [later] items stand.

Not pushed; no CI.

Taken in the recording commit: blocking item 3 (the status line).

## Ola's rulings, 2026-10-05 afternoon

Ola, verbatim: "yes to all defaults, push both". So: `tools/brief.py` takes the two edits (note names `<persona>-<worktree>-<HHMMSS>.md`, root `*.md` files in `@architect`'s limit); `count_loc.py` pins `--diff-algorithm=myers`; a file renamed across the counted/uncounted boundary counts in full (added on the way in, removed on the way out); `CLAUDE.md` §2's pointer says that `docs/` and non-code files are not counted either. Red for the second and third: `f885887`.

### Round 4: `@reviewer`, code, PR A whole, `bc01cd8..095ab63`

CHANGES REQUESTED. 287 net production lines (`count_loc.py bc01cd8 095ab63`). Round 3's items closed; 219 harness tests pass; rule lines match the tools. Blocking: the docstring of `test_env_drops_config_given_in_the_environment` names `diff.algorithm` as unpinned, which `195b504` made false; no CI before the push. Suggestion: cut "At 98e31cd only the one added line counts" from the rename test's docstring.

### Round 5: `@reviewer`, code, PR A, `095ab63..a6b966e`

CHANGES REQUESTED, on the status line only (taken in the recording commit). 287 net production lines (`count_loc.py bc01cd8 a6b966e`). Round 4's docstring items closed. Green CI after the push makes it APPROVED with no further round.

### Round 6: `@reviewer`, code, PR B, round 1

CHANGES REQUESTED. Recorded from the main session's summary in the brief
for this recording, not from the handback word for word: the handback was
in the previous session and was not carried over, and the summary does not
state the range (taken here to end at `084a6b3`, PR B's head when Ola
ruled) or a LOC count (78 net production lines by `tools/count_loc.py
origin/master 084a6b3`, measured for this record). Findings:

1. G3b (local git writes pass in a scratchpad repository): abbreviated or
   glued git options get through. `git -C <scratch repo> config --glo …`
   writes the global gitconfig, and `-f/path/to/real/.git/config` writes
   the real repository's config; the exemption matches exact spellings only.
2. G3a and G3b: re-pointing on the same command line. The guard judges the
   target before the line runs, so `ln -sfn <real repo> $SP/x && echo … >
   $SP/x/CLAUDE.md` passes, as does swapping a scratch repository's `.git`
   for a link first. A link made in an earlier command is already caught.
   [`$SP` stands for the scratchpad path written out: with a literal `$SP`
   the guard judges the static part, `/x/CLAUDE.md`, and asks.]
3. G1's fetch rule: `--refmap`, `--stdin` and `-c remote.<x>.fetch=…` write
   named refs with no `src:dst` on the command line.
4. Pre-existing: `gh -R owner/repo pr merge N` slips past the merge guard.

## Ola's ruling on review round 6, 2026-10-05

The options put to him: A, drop both exemptions (G3a and G3b); B, harden
both; C, keep G3a but only for a single plain command (no `&&`, `;`, `|`,
subshell and the like), and drop G3b entirely. The default offered with
them: fix the `gh -R` merge hole in PR B.

Ola, verbatim: "1C". So: option C, with the default. §2 G1 (amendment), G3
and G6, §3, §4, §5 and §6 are amended to match; finding 3 is fixed under
G1 whichever option was chosen.

### Round 7: `@reviewer`, code, PR B, round 2, `084a6b3..a42d864`

CHANGES REQUESTED. Range `084a6b3..a42d864` (design amendment 8554a3e, red 94989dd, green a42d864); 77 net production lines for PR B (`tools/count_loc.py origin/master a42d864`; merge base 7dde17a), against about 70. No CI yet (not pushed). Round 6's four findings closed, each probed through the hooks. Blocking: (1) `git -c include.path=<file> fetch s` (and `includeIf.*.path`, `--config-env=include.path=`) writes named refs with no question, probed with git 2.55; G1 asks only on `remote.`/`url.` keys, and §2 G1's "any other key … the repository's own config, which is governed" is false; (2) since h3, a glued `gh api -fquery=…`/`-Fquery=…` and curl `-d@f`/`-sd x`/`-Tf` to the forge pass, so a GraphQL `mergePullRequest` bypasses the merge guard; (3) the status line. Later: `find -exec` and `-c core.fsmonitor`/`core.hooksPath` named in §6; §7 question 3's default marked superseded. Pins accepted: bare `gh release`; `gh api` read from the first `api`.

Recorded word for word from `@reviewer`'s record text. Taken in the recording commit: blocking item 3 (the status line); blocking items 1 and 2 are §2 G1's and G6's round 7 amendments; the two later items are §6 and §7 question 3.

## Ola's ruling on review round 7, 2026-10-05

Ola, verbatim: "yes, close the gh api hole in PR B". So: PR B also closes the glued-option route past the merge and forge guard, `gh api … -fquery=…`/`-Fquery=…` and `curl -d@file`, `-sd x`, `-Tfile` to the forge (§2 G6, amendment after review round 7).

## Ola's ruling on §7 question 4, 2026-10-05

Ola, verbatim: "yes close it now, and defaults on the cell ID". The first
half answers question 4 (command runners); the second is about another
matter and is not ruled on here. So: PR B also closes the route past the
push guard through a program that runs its arguments as a command (§2 G7).
The governance guard's half stays in §6, for the reasons at the end of G7.

### Round 8: `@reviewer`, code, PR B, round 3, `a42d864..fa9f3e1`

Round 8: `@reviewer`, code, PR B, round 3, `a42d864..fa9f3e1`. CHANGES REQUESTED. The range is design 28b8292 + 1a1e5d3, red 8c309a8, green fa9f3e1. PR B is 97 net production lines (`tools/count_loc.py origin/master fa9f3e1`, merge base 7dde17a), against about 103; the round adds 20. Not pushed, so no CI. Locally: 768 passed in the four guard suites; ruff, ruff format, mypy, the prohibited-dependency and detria gates are green; the at-risk citations are all history and still read correctly; no red-step scaffolding is left. Round 7's items 1 and 2 are closed, each probed through the hook. Blocking:
- (1) The status line still says the red and green steps come next.
- (2) G7 rule (b) knows only `sh`, `bash` and `zsh` (`/Users/skavhaug/projects/rasputin/.claude/worktrees/h16b/tools/shell_scan.py:65` @fa9f3e1, copied at `/Users/skavhaug/projects/rasputin/.claude/worktrees/h16b/.claude/hooks/guard_push.py:87` @fa9f3e1). So `dash -c 'git push'`, `ksh`/`csh`/`tcsh -c …` and `caffeinate dash -c 'git push'` pass, and all four shells are in /bin on macOS.
- (3) `parallel git ::: push` passes: rule (a) reads the tail as `git ::: push`, and a bare tail drops the unknown-command reason (`/Users/skavhaug/projects/rasputin/.claude/worktrees/h16b/.claude/hooks/guard_push.py:307` @fa9f3e1).
- (4) A copy or move into a directory named without a trailing `/` is judged as the directory name (`/Users/skavhaug/projects/rasputin/.claude/worktrees/h16b/tools/shell_scan.py:359` @fa9f3e1, unchanged since h4). So `cp json.py tools`, `mv ast.py tools`, `cp s .git/remotes` and `cp x.py .claude/hooks` pass. PR B's new G4 rule and its new `.git/remotes/` prefix are both reachable this way.

Recorded word for word from `@reviewer`'s record text, taken from its handback in the session log, but for one path: item (3)'s short path to `guard_push.py`, line 307, is written out in full, as the record's other paths are, so `tools/check_citations.py` resolves it. Taken in the recording commit: blocking item 1 (the status line). Blocking items 2 and 3 are §2 G7's round 8 amendments, and item 4 is §2 G8. Later items: the `\;` strings in the tests go to the red step (§2 G7, *Red test, review round 8*), since the green commit touches no test file; the pre-existing routes the round listed are in §6.

## Ola's ruling on review round 8, 2026-10-06

The question, as the main session put it in its morning summary (row H, 2026-10-06 04:31 UTC): "Guard fix: close the 3 extra ways around the guards that the review found, in this PR (needs you here), or list them as known gaps?", with the default "close now".

Ola, verbatim: "P1 yes, all defaults" (04:34 UTC), and later "let's do the guard fix now" (11:06 UTC). So: PR B also closes the three routes review round 8 found: a git or gh command run by any shell (`dash`, `ksh`, `csh`, `tcsh` and the rest), a git or gh word under `parallel`, and a copy or move into a governed directory named without a trailing `/` (§2 G7 amended, G8 added).

### Round 9: `@reviewer`, code, PR B, round 4, `fe9e4511..a0001d00`

Round 9: `@reviewer`, code, PR B, round 4, `fe9e4511..a0001d00`. CHANGES REQUESTED. The range is the master merge f86d3953 (clean: `git show --remerge-diff` is empty), red 2d44b462 and green a0001d00. PR B is 109 net production lines (`tools/count_loc.py origin/master a0001d00`, merge base 879ea493), against about 106; the round adds 12 (`tools/count_loc.py f86d3953 a0001d00`). Not pushed, so no CI. Locally: 6125 passed and 28 skipped in the full Python suite, run on a scratch copy with the main checkout's venv; 1122 passed in the git-history suites and the guard suites in the worktree, with SyntaxWarning made an error; ruff, ruff format, mypy, the prohibited-dependency gate and the detria gate are green; no red-step scaffolding is left. Probed through the hooks, all of round 8's red rows ask and their controls and pins behave as designed. Blocking:
- (1) A runner in front of `parallel`, `watch` or `flock` hides it from rule (a)'s parallel exception and from rule (c). Both look only at the line's first program (`/Users/skavhaug/projects/rasputin/.claude/worktrees/h16b/.claude/hooks/guard_push.py:235` @a0001d00). So these pass: `caffeinate parallel git ::: push`, `stdbuf -o0 parallel git ::: push`, `arch -arm64 parallel git ::: push`, `find . -maxdepth 0 -exec parallel git ::: push \;`, `caffeinate watch 'git push'`, `caffeinate parallel 'git push' ::: a`, and `stdbuf -o0 flock /tmp/l -c 'git push'`. R1's *The harness* paragraph says a git or gh word under `parallel` asks.
- (2) G4's `tools/<standard-library name>` rule reads the path as written (`/Users/skavhaug/projects/rasputin/.claude/worktrees/h16b/.claude/hooks/guard_governance.py:116` @a0001d00). So `cp json.py tools/.`, `cp json.py tools/./` and `mv json.py tools/./json.py` pass. R1 says a copy into their directories asks.
- (3) The status line still says the merge, red and green steps come next.

Recorded word for word from `@reviewer`'s "One-line record" in its handback, taken from the session log. Taken in the recording commit: blocking item 3 (the status line). Blocking item 1 is §2 G7's round 9 amendment, and item 2 is G4's. Later items: the clustered `-t`, `--target=` and the whole-directory writes are in §6; the two stale line citations in `docs/increments/python-audit.md` (its table rows for this file's lines 579 and 615, which have moved) are pinned to `abd7c68`, where those lines read as cited.

## Ola's ruling on review round 9, 2026-10-06

The option the main session offered as its default, A: close the two routes now, and list writes to a whole governed directory (`rm -rf .claude/hooks`, `mv .claude/agents /tmp/`, `cp -R hooks .claude/`; also open on `master`) as a known gap for a later harness increment.

Ola, verbatim: "A, close the two now". So: PR B also closes a runner in front of `parallel`, `watch` or `flock` (§2 G7, amendment after review round 9) and a governed path written with `.`, `..` or a doubled `/` (§2 G4, amendment after review round 9). Writes to a whole governed directory are in §6. The slow line found while designing the G7 amendment is §7 question 5, not ruled on.

## Ola's ruling on §7 question 5, 2026-10-06

Ola, verbatim, 2026-10-06 22:13 UTC: "defaults on all four". So: the slow line is fixed in this round. The answer came before round 9's red step (`c83dfe0e`, 22:24 UTC), and round 9's green step (`80800a14`) fixed it.

### Round 10: `@reviewer`, code, PR B, round 5, `a0001d00..80800a14`

Round 10: `@reviewer`, code, PR B, round 5, `a0001d00..80800a14`. CHANGES REQUESTED. The range is red c83dfe0e and green 80800a14. PR B is 116 net production lines (`tools/count_loc.py origin/master HEAD`, merge base 879ea493) against about 115; the round adds 7 (`tools/count_loc.py a0001d00 HEAD`: guard_governance.py +4, guard_push.py +3). Not pushed, so no CI. Blockers: (1) the "linear" claim is false: every shell word parses the rest of the line again (.claude/hooks/guard_push.py@80800a14:238, :246; docs/increments/h16-harness-fixes.md@80800a14:764, :768), so time grows with shell words x line length (800 x sh plus 100 KB: 31 s); deny starts at about 990 sh words, not "a thousand or more"; (2) the status line (docs/increments/h16-harness-fixes.md@80800a14:3) and §7 question 5 (:1373) still read as unanswered. 1 suggestion: S1 an explicit cap on runner/shell words per line.

Recorded word for word from `@reviewer`'s record text, as the main session passed it on. Taken in the recording commit: blocking item 1 in this file (§2 G7's time and its constant); its other half, the docstring at `.claude/hooks/guard_push.py` line 238, is `@developer`'s, in its own commit. Blocking item 2 (the status line and §7 question 5). Suggestion S1 is in §6.

## Ola's ruling on review round 10, 2026-10-07

On suggestion S1, a cap on runner and shell words per line, the main session's default was no cap in this PR and a later harness increment. Ola, verbatim: "Ok, go for defaults." So: no cap now; §6 lists it, with the reason (it bounds the time, and replaces the crash-worded deny with a plain reason).

### Round 11: `@reviewer`, code, PR B, round 6, `80800a14..85d39f5c`

Round 11: `@reviewer`, code, PR B, round 6, `80800a14..85d39f5c`. APPROVED. The range is docs 19d6b3d7 and docstring 85d39f5c. PR B is still 116 net production lines (`tools/count_loc.py origin/master HEAD`), and the round adds 0 (`tools/count_loc.py 80800a14 HEAD`: guard_push.py 0). Not pushed, so no CI. Both round 10 blockers are fixed: (1) "linear" is gone from .claude/hooks/guard_push.py@85d39f5c:238 and docs/increments/h16-harness-fixes.md@85d39f5c:765, and the deny threshold is now about 990 `sh` words (:774). (2) The Status line (:3) and §7 question 5 (:1402) now record Ola's answer. 3 suggestions (S1 to S3).

Recorded word for word from `@reviewer`'s record line, as the main session passed it on. Its suggestions: S1, the `watch` count at §2 G7's constant (`docs/increments/h16-harness-fixes.md@85d39f5c:774-775`) holds only with the command after `watch` quoted; S2, the heading of Ola's rulings after round 10 held his answer to §7 question 5, which came before round 9's red step; S3, no change (the 31 s time is pinned to `80800a14`). Taken in the recording commit: S1 (§2 G7 now names each measured form, re-measured at `85d39f5c`: 993 quoted `watch` words here against the review's 995, 999 unquoted) and S2 (the ruling on question 5 has its own section, dated 2026-10-06, before round 10).
