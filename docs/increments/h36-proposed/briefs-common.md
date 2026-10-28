You are @$persona. This block comes from tools/brief.py; the task after it
is the main session's wording.
Read **.claude/agents/$persona.md from disk** and the increment file,
$increment. Where the task contradicts either, **the files win**; say so,
since **a brief cannot drop a step** they require.
Work only in $worktree (cd there first), with your own build directory and
venv; the write limit below is for it, and a scratch copy your persona file
names is yours to change and remove afterwards. Your note file is $note; write **no other file under .claude/current-task/**.
**Ola's words appear only under "Ola, verbatim"** below.
If **blocked on power, network or a lock**, stop and hand back.
Run a long job with Bash's `run_in_background` and wait for its end notice,
never `sleep` or a polling loop; before you hand back, stop each background
command you started, by its id.
Write a rule file with Edit or Write, never through the shell.
End each commit message with the **Co-Authored-By trailer** from your system
context, and each commit subject with `(@$persona)`.
Write in **plain words**: say what any internal label means.
Hand back, in at most 300 words, under: **Result; Pinned or assumed beyond the design; Questions
for Ola; Lessons; Ideas; DEFAULTS TAKEN and GUARD FALSE POSITIVE lines** ("none" under an
empty one). Each question for Ola is in plain words, with a default.
