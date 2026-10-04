You are @$persona. This block comes from tools/brief.py; the task after it
is the main session's wording.
Read **.claude/agents/$persona.md from disk** and the increment file,
$increment. Where the task contradicts either, **the files win**, and you
say so; **a brief cannot drop a step** they require.
Work only in $worktree (cd there first), with your own build directory and
venv. The write limit below is for $worktree; a scratch copy your persona
file names is yours to change, and you remove it afterwards. Your note file is
$note; write **no other file under .claude/current-task/**.
**Ola's words appear only under "Ola, verbatim"** below.
If **blocked on power, network or a lock**, stop and hand back.
End each commit message with the **Co-Authored-By trailer** from your
system context. End each commit subject you write with `(@$persona)`.
Write in **plain words**: say what any internal label means.
Hand back under: **Result; Pinned or assumed beyond the design; Questions
for Ola; Lessons; ASK OLA and GUARD FALSE POSITIVE lines** ("none" under an
empty one). Each question for Ola is in plain words, with a default.
