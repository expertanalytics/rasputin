"""What a shell line runs and which files it writes, for the guards (h4).

Spec: docs/increments/h4-guard-fixes.md §3. `parse(command)` splits a Bash line
into its simple commands, each with its argv (quotes removed), the files it
writes and, for an interpreter, the program it runs. Quotes, heredoc bodies and
comments are text, never commands, so a governed name or the words `git push`
inside them judge nothing. `$(...)`, backticks, `<(...)` and `sh -c S` are
parsed recursively; wrappers (`env`, `xargs`, `script -q F`, ...) are stripped.

Stdlib only and pure: the hooks run under whichever `python3` is on PATH. It is
a reader, not a shell: `None` means it could not read the line, and the guards
then judge the line as text, as before h4.
"""

from __future__ import annotations

import ast
import os
import re
from dataclasses import dataclass, field
from itertools import pairwise


@dataclass
class Simple:
    """One simple command, after wrappers and leading assignments are stripped."""

    argv: list[str]
    writes: list[str] = field(default_factory=list)
    #: The -c/-e text or the heredoc an interpreter reads; None for other commands.
    program: str | None = None
    #: The script (or `-m` module, as a path) an interpreter or shell runs.
    script: str | None = None
    #: It writes files it does not name (`git apply`, `patch`, `xargs rm`).
    unknown: bool = False
    #: The body of its heredoc, whoever reads it; `program` only for an interpreter.
    heredoc: str | None = field(default=None, repr=False)


class UnreadableError(Exception):
    """The line has unbalanced quotes, an unclosed `$(`, or an unterminated heredoc."""


RESERVED = {"if", "then", "else", "elif", "fi", "do", "done", "while", "until", "!", "{", "}"}
LOOP_HEADERS = {"for", "select", "function"}
ASSIGNMENT = re.compile(r"[A-Za-z_]\w*=")
REDIRECT = re.compile(r"(\d*|&)(>>|>\||>&|>|<<<|<<-|<<|<>|<&|<)")
WORD_END = set(" \t\n;&|()<>")
WRAPPERS = {"env", "nohup", "time", "timeout", "nice", "command", "exec", "xargs", "script", "sudo"}
#: Options of the wrappers above whose next word is their argument.
WRAPPER_ARGS = {"-u", "-n", "-s", "-k", "-g", "-I", "-P", "-L", "-d", "-E"}
#: Writers by argv, with the options whose next word is their argument.
WRITER_ARGS = {
    "tee": set(),
    "rm": set(),
    "rmdir": set(),
    "unlink": set(),
    "mv": {"-t", "-S"},
    "touch": {"-t", "-r", "-d"},
    "truncate": {"-s", "-r"},
    "mkdir": {"-m"},
}
COPIERS = {"cp": {"-t", "-S"}, "ln": {"-t", "-S"}, "install": {"-m", "-o", "-g", "-t", "-S"}}
COPIERS["mv"] = WRITER_ARGS.pop("mv")  # h16 G8: read as a copier that also names its sources
PYTHON = re.compile(r"python(3(\.\d+)?)?")
#: Shells that take a program with -c: /etc/shells on macOS, plus common ones elsewhere (h16 G7).
SHELLS = {"sh", "bash", "zsh", "dash", "ksh", "csh", "tcsh", "rbash", "fish", "mksh", "ash", "yash"}
#: Interpreter -> the option letter that introduces its program text.
PROGRAM_FLAG = {"perl": "e", "ruby": "e", "node": "e"} | dict.fromkeys(SHELLS, "c")
PLAIN_PARAMETER = re.compile(r"\$(?:\{(\w+)\}|(\w+))")
PARAMETER = re.compile(r"\$[\w@*#?$!-]?\w*")
EXPANSION = re.compile(r"\$(\{[^}]*\}|\([^)]*\)|[\w@*#?$!-]+)|`[^`]*`")


def base(word: str) -> str:
    return word.rsplit("/", 1)[-1]


def static(target: str, filler: str = "") -> str:
    """The target with every expansion replaced by `filler`: `$D/CLAUDE.md` -> `/CLAUDE.md`."""
    return EXPANSION.sub(filler, target)


def ignored(target: str) -> bool:
    return target == "/dev/null" or target.startswith(("/dev/std", "/dev/fd/"))


def parse(command: str) -> list[Simple] | None:
    """The simple commands of `command`, in the order the reader meets them, or None."""
    scan = _Scan(command)
    try:
        scan.run(None)
    except UnreadableError:
        return None
    return scan.found


class _Scan:
    def __init__(self, text: str) -> None:
        self.text, self.at = text, 0
        self.found: list[Simple] = []
        self.pending: list[tuple[Simple, str, bool]] = []  # heredocs: owner, delimiter, <<-
        self.backtick = False
        self.assigned: dict[str, str] = {}  # NAME=value commands met so far on the line

    def run(self, close: str | None) -> None:
        """Read commands up to `close` (`)` or a backtick) or, at the top, the end."""
        words: list[str] = []
        writes: list[str] = []
        heredocs: list[tuple[str, bool]] = []
        while True:
            self.blanks()
            at, text = self.at, self.text
            if at >= len(text) or (close is not None and text[at] == close):
                if close is not None and at >= len(text):
                    raise UnreadableError
                self.at += close is not None
                self.finish(words, writes, heredocs)
                if close is None and self.pending:
                    raise UnreadableError
                return
            char = text[at]
            if REDIRECT.match(text, at):
                if char in "<>" and text[at + 1 : at + 2] == "(":
                    words.append(self.word())  # process substitution
                    continue
                self.redirect(writes, heredocs)
            elif char in "\n;&|":
                self.at += 2 if text[at : at + 2] in ("&&", "||", ";;", "|&") else 1
                self.finish(words, writes, heredocs)
                words, writes, heredocs = [], [], []
                if char == "\n":
                    self.read_heredocs()
            elif char == "(":
                self.at += 1
                if words and text[self.at : self.at + 1] == ")":
                    self.at += 1  # a function definition, `f()`
                    continue
                self.finish(words, writes, heredocs)
                words, writes, heredocs = [], [], []
                self.run(")")
            elif char == ")":
                raise UnreadableError
            else:
                words.append(self.word())

    def blanks(self) -> None:
        text = self.text
        while self.at < len(text):
            if text[self.at] in " \t":
                self.at += 1
            elif text.startswith("\\\n", self.at):
                self.at += 2
            elif text[self.at] == "#":
                end = text.find("\n", self.at)
                self.at = len(text) if end < 0 else end
            else:
                return

    def word(self) -> str:
        text, start = self.text, self.at
        out: list[str] = []
        if text[self.at] in "<>":  # <(...) or >(...)
            self.at += 2
            self.run(")")
            return text[start : self.at]
        while self.at < len(text):
            char = text[self.at]
            if char in WORD_END or (char == "`" and self.backtick):
                break
            if char == "'":
                end = text.find("'", self.at + 1)
                if end < 0:
                    raise UnreadableError
                out.append(text[self.at + 1 : end])
                self.at = end + 1
            elif char == '"':
                out.append(self.quoted())
            elif char == "\\":
                out.append(text[self.at + 1 : self.at + 2].replace("\n", ""))
                self.at += 2
            else:
                out.append(self.dollar() if char in "$`" else char)
                self.at += char not in "$`"
        if self.at == start:
            raise UnreadableError  # a backtick where none can close
        return "".join(out)

    def quoted(self) -> str:
        """A double-quoted string from its opening quote; substitutions are parsed."""
        text, out = self.text, list[str]()
        self.at += 1
        while self.at < len(text):
            char = text[self.at]
            if char == '"':
                self.at += 1
                return "".join(out)
            if char == "\\" and text[self.at + 1 : self.at + 2] in ('"', "\\", "$", "`", "\n"):
                out.append(text[self.at + 1].replace("\n", ""))
                self.at += 2
            elif char in "$`":
                out.append(self.dollar())
            else:
                out.append(char)
                self.at += 1
        raise UnreadableError

    def dollar(self) -> str:
        """A `$` or backtick construct, as written; its commands are added to `found`."""
        text, start = self.text, self.at
        if text[start] == "`":
            self.at, outer, self.backtick = start + 1, self.backtick, True
            self.run("`")
            self.backtick = outer
        elif text.startswith("$(", start):
            self.at = start + 2
            self.run(")")
        elif text.startswith(("${", "$'"), start):
            close = "}" if text[start + 1] == "{" else "'"
            end = text.find(close, start + 2)
            if end < 0:
                raise UnreadableError
            self.at = end + 1
        else:
            match = PARAMETER.match(text, start)
            assert match is not None
            self.at = match.end()
        return text[start : self.at]

    def redirect(self, writes: list[str], heredocs: list[tuple[str, bool]]) -> None:
        match = REDIRECT.match(self.text, self.at)
        assert match is not None
        self.at = match.end()
        prefix, op = match.groups()
        self.blanks()
        if self.at >= len(self.text) or self.text[self.at] in WORD_END:
            raise UnreadableError
        target = self.word()
        if op in ("<<", "<<-"):  # quoted or not: a body's expansions are not read here
            heredocs.append((target, op == "<<-"))
        elif op == ">&" and (target.isdigit() or target == "-"):
            return
        elif (op in (">", ">>", ">|", "<>", ">&") or prefix == "&") and not ignored(target):
            writes.append(target)

    def read_heredocs(self) -> None:
        """Consume the bodies of the heredocs opened on the line just ended."""
        text = self.text
        for owner, delimiter, tabs in self.pending:
            body: list[str] = []
            while True:
                if self.at >= len(text):
                    raise UnreadableError
                end = text.find("\n", self.at)
                end = len(text) if end < 0 else end
                line = text[self.at : end]
                self.at = min(end + 1, len(text))
                if (line.lstrip("\t") if tabs else line) == delimiter:
                    break
                body.append(line + "\n")
            owner.heredoc = "".join(body)
            self.finish_heredoc(owner)
        self.pending = []

    def finish(self, words: list[str], writes: list[str], heredocs: list[tuple[str, bool]]) -> None:
        if words and all(ASSIGNMENT.match(word) for word in words):  # `H=.git/harness;`
            self.assigned |= dict(word.split("=", 1) for word in words)
        while words and (words[0] in RESERVED or ASSIGNMENT.match(words[0])):
            words = words[1:]
        if words and words[0] in LOOP_HEADERS:
            words = []
        if not words and not writes and not heredocs:
            return
        simple = classify(words, writes)
        # A variable assigned earlier on the line is read as its value: `$H/x`.
        simple.writes = [PLAIN_PARAMETER.sub(self.value, w) for w in simple.writes]
        self.pending += [(simple, delimiter, tabs) for delimiter, tabs in heredocs]
        if simple.argv[:1] == ["eval"]:  # eval runs its words as a command line
            self.nested(" ".join(simple.argv[1:]))
            if not simple.writes:
                return
        self.found.append(simple)
        if simple.argv and base(simple.argv[0]) in SHELLS and simple.program is not None:
            self.nested(simple.program)
            simple.program = None

    def value(self, match: re.Match[str]) -> str:
        return self.assigned.get(match.group(1) or match.group(2), match.group(0))

    def finish_heredoc(self, simple: Simple) -> None:
        """A heredoc is the program of an interpreter reading stdin; a shell's is parsed."""
        if not simple.argv or simple.script not in (None, "-") or simple.program is not None:
            return
        name = base(simple.argv[0])
        if name in SHELLS:
            self.nested(simple.heredoc or "")
        elif PYTHON.fullmatch(name) or name in PROGRAM_FLAG:
            simple.program = simple.heredoc

    def nested(self, command: str) -> None:
        inner = parse(command)
        if inner is None:
            raise UnreadableError
        self.found += inner


def classify(words: list[str], writes: list[str]) -> Simple:
    """Strip wrappers, then name what the command writes and, if it interprets, its program."""
    argv, wrapped, stdin_args = unwrap(words)
    simple = Simple(argv, [w for w in [*writes, *wrapped] if not ignored(w)])
    if not argv:
        return simple
    name, args = base(argv[0]), argv[1:]
    named, unknown = writer_targets(name, args)
    simple.writes += named
    simple.unknown = unknown or (stdin_args and name in {*WRITER_ARGS, *COPIERS, "sed"})
    if PYTHON.fullmatch(name) or name in PROGRAM_FLAG:
        simple.program, simple.script = interpreter(name, args)
    return simple


def unwrap(argv: list[str]) -> tuple[list[str], list[str], bool]:
    """(argv without wrappers, files the wrappers write, whether xargs supplies arguments)."""
    writes, stdin_args = [], False
    while argv and base(argv[0]) in WRAPPERS:
        name, rest = base(argv[0]), argv[1:]
        while rest and rest[0].startswith("-") and rest[0] != "-":
            rest = rest[2:] if rest[0] in WRAPPER_ARGS else rest[1:]
        while name == "env" and rest and ASSIGNMENT.match(rest[0]):
            rest = rest[1:]
        if name in ("timeout", "script") and rest:
            writes += rest[:1] if name == "script" else []  # script's typescript file
            rest = rest[1:]
        stdin_args |= name == "xargs"
        argv = rest
    return argv, writes, stdin_args


def operands(args: list[str], takes: set[str]) -> list[str]:
    found: list[str] = []
    skip = False
    for at, arg in enumerate(args):
        if skip:
            skip = False
        elif arg == "--":
            return found + args[at + 1 :]
        elif arg in takes:
            skip = True
        elif not arg.startswith("-") or arg == "-":
            found.append(arg)
    return found


def writer_targets(name: str, args: list[str]) -> tuple[list[str], bool]:
    """(the files a writer names, whether it also writes files it does not name)."""
    if name in WRITER_ARGS:
        return operands(args, WRITER_ARGS[name]), False
    if name in COPIERS:  # h16 G8: a target directory, by option or by a name with no extension
        args = [word for arg in args for word in target_option(arg)]
        found, into = operands(args, COPIERS[name]), [b for a, b in pairwise(args) if a == "-t"]
        kept = found if name == "mv" else []  # mv also removes its sources
        if into:
            return kept + [f"{into[-1].rstrip('/')}/{base(source)}" for source in found], False
        if len(found) < 2:
            return found, False
        dest = found[-1]
        inside = [f"{dest.rstrip('/')}/{base(source)}" for source in found[:-1]]
        if dest in (".", "..") or dest.endswith("/"):  # into a directory: the sources' names
            return kept[:-1] + inside, False
        return [*kept[:-1], dest, *(inside if not os.path.splitext(dest)[1] else [])], False
    if name == "sed" and any(a.startswith(("-i", "--in-place")) for a in args):
        return sed_files(args), False
    if name == "dd":
        return [a[3:] for a in args if a.startswith("of=")], False
    if name == "patch":
        return [], True
    if name == "git":
        at = 0
        while at < len(args) and args[at].startswith("-"):
            at += 2 if args[at] in ("-C", "-c") else 1
        sub, rest = (args[at], args[at + 1 :]) if at < len(args) else ("", [])
        if sub in ("checkout", "restore", "mv", "rm"):
            return operands(rest, {"-b", "-B", "--source", "-s"}), False
        return [], sub in ("apply", "am")
    return [], False


def target_option(arg: str) -> list[str]:
    """GNU's `-tDIR`, `--target-directory DIR` and `--target-directory=DIR`, as `-t DIR`."""
    if arg.startswith("--target-directory"):
        return ["-t", *arg.split("=", 1)[1:]]
    return ["-t", arg[2:]] if arg.startswith("-t") and len(arg) > 2 else [arg]


def sed_files(args: list[str]) -> list[str]:
    """The files `sed -i` edits: every operand but the script, `-i ''` included."""
    found: list[str] = []
    script_given, at = False, 0
    while at < len(args):
        arg = args[at]
        if arg == "-i" and args[at + 1 : at + 2] == [""]:
            at += 1  # BSD sed's empty backup suffix
        elif arg in ("-e", "-f"):
            script_given, at = True, at + 1
        elif not arg.startswith("-"):
            found.append(arg)
        at += 1
    return found if script_given else found[1:]


def interpreter(name: str, args: list[str]) -> tuple[str | None, str | None]:
    """(program text, script operand) of an interpreter; `-m mod` is the script `mod.py`."""
    letter = "c" if PYTHON.fullmatch(name) else PROGRAM_FLAG[name]
    at = 0
    while at < len(args):
        arg = args[at]
        if arg == "-" or not arg.startswith("-"):
            return None, arg
        if arg.startswith("--"):
            if arg == "--eval" and at + 1 < len(args):
                return args[at + 1], None
        elif arg.endswith(letter) and at + 1 < len(args):
            return args[at + 1], None
        elif arg == "-m" and PYTHON.fullmatch(name) and at + 1 < len(args):
            return None, args[at + 1].replace(".", "/") + ".py"
        elif arg in ("-W", "-X"):
            at += 1
        at += 1
    return None, None


def candidates(simple: Simple) -> list[str]:
    """Files an interpreter command may write, read from its program and operands.

    Python: the string literals that are a single path (no whitespace); prose
    literals are not paths. Other languages, or Python that does not parse: the
    program's path-like tokens, and for those other languages the operands too
    (`perl -pi -e S FILE`).
    """
    name = base(simple.argv[0]) if simple.argv else ""
    if not (PYTHON.fullmatch(name) or name in PROGRAM_FLAG) or name in SHELLS:
        return []
    program = simple.program or ""
    if PYTHON.fullmatch(name):
        try:
            tree = ast.parse(program)
        except (SyntaxError, ValueError):
            return re.findall(r"[^\s'\"(),;<>{}\[\]=]+", program)
        return [
            node.value
            for node in ast.walk(tree)
            if isinstance(node, ast.Constant)
            and isinstance(node.value, str)
            and node.value
            and not any(c.isspace() for c in node.value)
        ]
    tokens = re.findall(r"[^\s'\"(),;<>{}\[\]=]+", program)
    return tokens + operands(simple.argv[1:], set())
