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

A command that could run text (a runner) whose stdin or script could carry text
the reader did not read (it is fed) is `unread`, and so `unknown`: the guards
judge its line by its text. A tripwire, not a sandbox; what stays open is in
`docs/increments/h25-stdin-routes.md`, *3.6 What stays uncovered*.
"""

from __future__ import annotations

import ast
import builtins
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
    #: A runner fed text the reader did not read; it is then `unknown` too.
    unread: bool = False
    #: Its words before wrappers are stripped; whether it stands after a pipe; whether
    #: its own command has an input redirect or a process substitution.
    written: list[str] = field(default_factory=list, repr=False, compare=False)
    piped: bool = field(default=False, repr=False, compare=False)
    fed: bool = field(default=False, repr=False, compare=False)


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
#: The names a read relies on, which an exempt program may not rebind.
READ_NAMES = {"open", "json", "pathlib", "Path"}
#: Names that patch, reflect or run code; an exempt program uses none of them.
REFLECTION = {"setattr", "delattr", "getattr", "vars", "locals", "globals", "exec", "eval",
              "compile", "breakpoint", "__import__", "__builtins__", "builtins"}  # fmt: skip
EXPANSION = re.compile(r"\$(\{[^}]*\}|\([^)]*\)|[\w@*#?$!-]+)|`[^`]*`")
#: Command words that run text: their own, and among any command's words (Python by PYTHON).
RUNNER_FIRST = {".", "source", "eval", "exec", "coproc", "parallel"}
RUNNER_ANY = {*SHELLS, "perl", "ruby", "node", "uv", "uvx"}
#: A word token naming a runner, bare or after an option, `=` or `/` (`-Spython3`, `/bin/sh`).
RUNNER_TOKEN = re.compile(rf"(?:-\w*|\S*[=/])?(?:{'|'.join(RUNNER_ANY)}|{PYTHON.pattern})")
#: Commands that run no program whatever their arguments: a runner's name among them is data.
NO_RUN = {"grep", "egrep", "fgrep", "rg", "head", "tail", "cat", "wc", "sort", "uniq", "cut", "tr",
          "ls", "echo", "printf", "diff", "cmp", "comm", "basename", "dirname", "realpath", "test",
          "[", "which", "type", "file", "stat", "du", "df", "tee", "mkdir", "rm", "cp", "mv",
          "touch"}  # fmt: skip
OPENERS = {"{", "if", "while", "until", "for", "case", "select"}
CLOSERS = {"}", "fi", "done", "esac"}
#: A command word that a glob, a brace or zsh's `=name` could turn into any program.
PATTERN = re.compile(r"[*?\[{]|^=")
#: The allowlist of a Python data reader: flags, modules and what a program may touch.
PLAIN_FLAGS = set("BEIOPqsSubdRvx")
PLAIN_PATH = re.compile(r"[^$`*?\[{<>=()]+")
DATA_MODULES = {"sys": {"stdin", "stdout", "stderr", "argv", "exit"},
                "json": {"load", "loads", "dump", "dumps"}}  # fmt: skip
DATA_BUILTINS = {"print", "len", "sorted", "min", "max", "sum", "any", "all", "abs", "round",
                 "int", "float", "str", "bool", "list", "dict", "set", "tuple", "enumerate", "zip",
                 "range", "map", "filter", "reversed", "isinstance", "repr", "iter", "next",
                 "divmod", "KeyError", "ValueError", "IndexError", "Exception"}  # fmt: skip
DATA_METHODS = {"read", "readline", "readlines", "write", "strip", "lstrip", "rstrip", "split",
                "rsplit", "splitlines", "join", "replace", "startswith", "endswith", "lower",
                "upper", "get", "items", "keys", "values", "append", "extend", "pop", "sort",
                "count", "index", "find", "partition", "rpartition", "setdefault", "update",
                "format", "ljust", "rjust", "isdigit"}  # fmt: skip
NOT_DATA = (ast.ImportFrom, ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef, ast.Global,
            ast.Nonlocal, ast.Await, ast.Yield, ast.YieldFrom)  # fmt: skip


def base(word: str) -> str:
    return word.rsplit("/", 1)[-1]


def static(target: str, filler: str = "") -> str:
    """The target with every expansion replaced by `filler`: `$D/CLAUDE.md` -> `/CLAUDE.md`."""
    return EXPANSION.sub(filler, target)


def ignored(target: str) -> bool:
    return target == "/dev/null" or target.startswith(("/dev/std", "/dev/fd/"))


def parse(command: str, piped: bool = False) -> list[Simple] | None:
    """The simple commands of `command`, in the order the reader meets them, or None.

    `piped`: the line reads a pipe, as a piped shell's `-c` program or a shell's heredoc does.
    """
    scan = _Scan(command)
    try:
        scan.run(None, piped)
    except UnreadableError:
        return None
    feeders = [s for s in scan.heredoc_owners if runner(s, scan.defined)]
    for simple in scan.found:
        fed = simple.piped or simple.fed or scan.line_fed or bool(feeders) or by_variable(simple)
        if fed and runner(simple, scan.defined) and not data_reader(simple):
            simple.unread = simple.unknown = True
    return scan.found


class _Scan:
    def __init__(self, text: str) -> None:
        self.text, self.at = text, 0
        self.found: list[Simple] = []
        # heredocs: owner, delimiter, <<-, an unquoted delimiter
        self.pending: list[tuple[Simple, str, bool, bool]] = []
        self.backtick = False
        self.assigned: dict[str, str] = {}  # NAME=value commands met so far on the line
        # Stdin state of the command being read: after a pipe; the piped state each open
        # compound restores at `;` and newlines; an input redirect; a pipe with nothing after.
        self.piped, self.stack, self.inputs, self.fresh = False, [False], False, False
        self.line_fed = False  # a feed on a compound, `exec` or `coproc` reaches every runner
        self.defined: set[str] = set()  # functions and aliases defined on the line
        self.heredoc_owners: list[Simple] = []  # bodies on a command not literally a reader

    def run(self, close: str | None, inherit: bool = False) -> None:
        """Read commands up to `close` (`)` or a backtick) or, at the top, the end."""
        outer = self.piped, self.stack, self.inputs, self.fresh
        self.piped, self.stack, self.inputs, self.fresh = inherit, [inherit], False, False
        words: list[str] = []
        writes: list[str] = []
        heredocs: list[tuple[str, bool | None, bool]] = []
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
                self.line_fed |= self.fresh  # zsh's trailing `|&` starts a coprocess
                self.piped, self.stack, self.inputs, self.fresh = outer
                return
            char = text[at]
            fresh, self.fresh = self.fresh, False
            if REDIRECT.match(text, at):
                if char in "<>" and text[at + 1 : at + 2] == "(":
                    words.append(self.word())  # process substitution
                    continue
                self.redirect(writes, heredocs)
            elif char in "\n;&|":
                sep = text[at : at + 2] if text[at : at + 2] in ("&&", "||", ";;", "|&") else char
                self.at += len(sep)
                self.finish(words, writes, heredocs)
                words, writes, heredocs = [], [], []
                if sep in ("|", "|&"):
                    self.piped = self.fresh = True
                elif char == "\n" and fresh:  # a newline right after a pipe continues it
                    self.fresh = True
                else:
                    self.piped = self.stack[-1]
                if char == "\n":
                    self.read_heredocs()
            elif char == "(":
                self.at += 1
                if words and text[self.at : self.at + 1] == ")":
                    self.at += 1  # a function definition, `f()`
                    self.defined.add(words[-1])
                    continue
                self.finish(words, writes, heredocs)
                words, writes, heredocs = [], [], []
                self.run(")", self.piped)
            elif char == ")":
                raise UnreadableError
            else:
                word, head = self.word(), all(w in RESERVED for w in words)
                if head and word in OPENERS:
                    self.stack = [*self.stack, self.piped]
                elif head and word in CLOSERS and len(self.stack) > 1:
                    self.stack = self.stack[:-1]
                words.append(word)

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
        if text[self.at] in "<>=" and text[self.at + 1 : self.at + 2] == "(":  # <(…) >(…) =(…)
            self.at += 2
            self.run(")", self.piped or text[start] == ">")
            self.inputs = True
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
            self.run("`", self.piped)
            self.backtick = outer
        elif text.startswith("$(", start):
            self.at = start + 2
            self.run(")", self.piped)
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

    def redirect(self, writes: list[str], heredocs: list[tuple[str, bool | None, bool]]) -> None:
        match = REDIRECT.match(self.text, self.at)
        assert match is not None
        self.at = match.end()
        prefix, op = match.groups()
        self.blanks()
        substitution = self.text.startswith(("<(", ">("), self.at)  # never a file written
        if self.at >= len(self.text) or (self.text[self.at] in WORD_END and not substitution):
            raise UnreadableError
        start, target = self.at, self.word()
        if op in ("<", "<&", "<>") and target not in ("/dev/null", "-"):
            self.inputs = True
        if op in ("<<", "<<-", "<<<"):  # a here-string (`<<<`) has None for `<<-`
            expands = not any(c in self.text[start : self.at] for c in "'\"\\")
            heredocs.append((target, None if op == "<<<" else op == "<<-", expands))
        elif substitution or (op == ">&" and (target.isdigit() or target == "-")):
            return
        elif (op in (">", ">>", ">|", "<>", ">&") or prefix == "&") and not ignored(target):
            writes.append(target)

    def read_heredocs(self) -> None:
        """Consume the bodies of the heredocs opened on the line just ended."""
        text = self.text
        for owner, delimiter, tabs, expands in self.pending:
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
            rewritten = expands and any(c in owner.heredoc for c in "$`\\")
            self.finish_heredoc(owner, owner.heredoc, rewritten)
        self.pending = []

    def finish(
        self, words: list[str], writes: list[str], heredocs: list[tuple[str, bool | None, bool]]
    ) -> None:
        inputs, self.inputs = self.inputs, False
        closed = not words or words[0] in CLOSERS  # a compound's or a subshell's redirects
        if words and all(ASSIGNMENT.match(word) for word in words):  # `H=.git/harness;`
            self.assigned |= dict(word.split("=", 1) for word in words)
        while words and (words[0] in RESERVED or ASSIGNMENT.match(words[0])):
            words = words[1:]
        if words[:1] == ["function"] and len(words) > 1:
            self.defined.add(words[1])
        if words and words[0] in LOOP_HEADERS:
            words = []
        head = words[0] if words else ""
        if ((closed or head == "exec") and (inputs or heredocs)) or head == "coproc":
            self.line_fed = True
        if not words and not writes and not heredocs:
            return
        simple = classify(words, writes)
        simple.written, simple.piped, simple.fed = words, self.piped, inputs
        if simple.argv[:1] == ["alias"]:
            self.defined |= {word.split("=", 1)[0] for word in simple.argv[1:] if "=" in word}
        # A variable assigned earlier on the line is read as its value: `$H/x`.
        simple.writes = [PLAIN_PARAMETER.sub(self.value, w) for w in simple.writes]
        self.pending += [(simple, d, tabs, x) for d, tabs, x in heredocs if tabs is not None]
        if simple.argv[:1] == ["eval"]:  # eval runs its words as a command line
            self.nested(" ".join(simple.argv[1:]), self.piped)
            if not simple.writes and not by_variable(simple):
                return
        self.found.append(simple)
        if simple.argv and base(simple.argv[0]) in SHELLS and simple.program is not None:
            self.nested(simple.program, self.piped)
            simple.program = None
        for word, tabs, _ in heredocs:
            if tabs is None:  # a here-string: the word and a newline, on stdin
                self.finish_heredoc(simple, word + "\n", substituted(word))

    def value(self, match: re.Match[str]) -> str:
        return self.assigned.get(match.group(1) or match.group(2), match.group(0))

    def finish_heredoc(self, simple: Simple, body: str, rewritten: bool) -> None:
        """A stdin body is parsed by a shell, and is the program of an interpreter reading stdin.

        Only for a command word that is literally an interpreter, a shell, `.` or `source`;
        on any other command the body feeds every runner on the line, if it is a runner.
        A shell's body is parsed as piped commands, whatever its `-c` program, since that
        program may read stdin as commands (`bash -c bash`); a shell with any argument, or
        a body that does not parse, feeds the line. An interpreter's body is appended
        whatever its options, as zsh joins them all, unless it is a data reader;
        `rewritten` marks a body the shell may rewrite first, which makes it unknown.
        """
        name = literal_name(simple)
        if not (PYTHON.fullmatch(name) or name in PROGRAM_FLAG or name in (".", "source")):
            self.heredoc_owners.append(simple)
        elif name in SHELLS or name in (".", "source"):
            inner = parse(body, True)
            self.line_fed |= inner is None or len(simple.argv) > 1
            self.found += inner or []
        elif not data_reader(simple):
            program = simple.program or ""
            simple.program = program + ("\n" if program and program[-1] != "\n" else "") + body
            simple.unknown = simple.unknown or rewritten

    def nested(self, command: str, piped: bool) -> None:
        inner = parse(command, piped)
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
        rewritten = simple.program is not None and substituted(simple.program)
        simple.unknown = simple.unknown or (name not in SHELLS and rewritten)
    return simple


def literal_name(simple: Simple) -> str:
    """The command word's name, after wrappers, when it holds no expansion or pattern; else ""."""
    first = simple.argv[0] if simple.argv else ""
    return "" if substituted(first) or PATTERN.search(first) else base(first)


def runner(simple: Simple, defined: set[str]) -> bool:
    """It could run text: a runner's name, a word that could be any program, or one defined."""
    written, argv = simple.written, simple.argv
    if not written:
        return False
    names = {base(written[0]), *(base(word) for word in argv[:1])}
    if names & (RUNNER_FIRST | RUNNER_ANY | defined) or any(PYTHON.fullmatch(n) for n in names):
        return True
    heads = {written[0], *argv[:1]} - {"[", "[[", "{", "}"}
    if any(substituted(word) or PATTERN.search(word) for word in heads):
        return True
    tokens = [token for word in written for token in word.split()]
    if not (argv and base(argv[0]) in NO_RUN) and any(map(RUNNER_TOKEN.fullmatch, tokens)):
        return True
    xargs = any(base(word) == "xargs" for word in written)
    replaces = ("{}" in w or w.startswith(("-I", "-J", "-i", "--replace")) for w in written)
    return xargs and (any(replaces) or any(map(substituted, written)))


def by_variable(simple: Simple) -> bool:
    """Fed through a variable: `eval`, `.`, `source` or an inline program holding an expansion."""
    name = base(simple.argv[0]) if simple.argv else ""
    if name in ("eval", ".", "source"):
        return any(map(substituted, simple.argv[1:]))
    if PYTHON.fullmatch(name) or name in PROGRAM_FLAG:
        program = interpreter(name, simple.argv[1:])[0]
        return program is not None and substituted(program)
    return False


def data_reader(simple: Simple) -> bool:
    """Plain `python3` with no-value flags and a plain script path or a data program."""
    words = simple.written
    if not words or not PYTHON.fullmatch(words[0]):
        return False
    at = 1
    while at < len(words) and words[at].startswith("-"):
        if words[at] == "-c":
            return at + 1 < len(words) and data_program(words[at + 1])
        if len(words[at]) < 2 or not set(words[at][1:]) <= PLAIN_FLAGS:
            return False
        at += 1
    script = words[at] if at < len(words) else ""
    if not PLAIN_PATH.fullmatch(script) or script.startswith("//"):
        return False
    path = os.path.normpath(script)  # `.` and `..` resolved: `/usr/../dev/stdin` is stdin
    return ".." not in path.split("/") and not path.startswith("/dev/")


def data_program(text: str) -> bool:
    """A Python program that only reads and writes data: every node is on the allowlist.

    It imports only `json` and `sys`, names only those, the builtins of DATA_BUILTINS
    and names it binds itself (never a builtin's, never one starting `_`), and reads
    only the attributes of DATA_MODULES and the methods of DATA_METHODS.
    """
    try:
        tree = ast.parse(text)
    except (SyntaxError, ValueError):
        return False
    nodes = list(ast.walk(tree))
    bound = {n.id for n in nodes if isinstance(n, ast.Name) and not isinstance(n.ctx, ast.Load)}
    bound |= {n.arg for n in nodes if isinstance(n, ast.arg)}
    bound |= {n.name for n in nodes if isinstance(n, ast.ExceptHandler) and n.name}
    if substituted(text) or bound & {*DATA_MODULES, *dir(builtins)}:
        return False
    allowed = {*DATA_MODULES, *DATA_BUILTINS, *bound}
    for node in nodes:
        if isinstance(node, NOT_DATA):
            return False
        names = node.names if isinstance(node, ast.Import) else []
        if any(alias.name not in DATA_MODULES or alias.asname for alias in names):
            return False
        if isinstance(node, ast.Name) and (node.id.startswith("_") or node.id not in allowed):
            return False
        if isinstance(node, ast.Attribute):
            owner = node.value.id if isinstance(node.value, ast.Name) else ""
            readable = DATA_MODULES.get(owner, DATA_METHODS)
            if not isinstance(node.ctx, ast.Load) or node.attr not in readable:
                return False
    return True


def substituted(text: str) -> bool:
    """The shell may splice text into it: it holds `$` or a backtick."""
    return "$" in text or "`" in text


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
        exempt = set() if simple.unknown else read_literals(tree)
        return [
            node.value
            for node in ast.walk(tree)
            if isinstance(node, ast.Constant)
            and isinstance(node.value, str)
            and node.value
            and (substituted(node.value) or not any(c.isspace() for c in node.value))
            and id(node) not in exempt
        ]
    tokens = re.findall(r"[^\s'\"(),;<>{}\[\]=]+", program)
    return tokens + operands(simple.argv[1:], set())


def literal(node: ast.expr, among: tuple[str, ...] = ()) -> bool:
    """A string literal; one of `among` when that is given."""
    found = isinstance(node, ast.Constant) and isinstance(node.value, str)
    return found and (not among or getattr(node, "value", None) in among)


def literal_keywords(call: ast.Call) -> bool:
    """Every keyword is `encoding` or `errors`, given as a string literal."""
    return all(k.arg in ("encoding", "errors") and literal(k.value) for k in call.keywords)


def opened(node: ast.expr, name: str = "open") -> ast.expr | None:
    """LIT of `open(LIT)` or `open(LIT, M)` with M a read mode (or `Path(LIT)`), else None."""
    if not isinstance(node, ast.Call) or not node.args or not literal(node.args[0]):
        return None
    func, modes = node.func, node.args[1:]
    if name == "open":
        ok = len(modes) < 2 and all(literal(m, ("r", "rb", "rt")) for m in modes)
        ok = ok and isinstance(func, ast.Name) and func.id == "open" and literal_keywords(node)
    else:
        named = isinstance(func, ast.Name) and func.id == "Path"
        dotted = isinstance(func, ast.Attribute) and func.attr == "Path"
        owner = getattr(func, "value", None)
        dotted = dotted and isinstance(owner, ast.Name) and owner.id == "pathlib"
        ok = not modes and not node.keywords and (named or dotted)
    return node.args[0] if ok else None


def read_of(node: ast.AST) -> ast.expr | None:
    """The literal a call reads and nothing else reaches: the three exempt shapes."""
    if not isinstance(node, ast.Call) or not isinstance(node.func, ast.Attribute):
        return None
    func, plain = node.func, not node.args and not node.keywords
    if func.attr == "load" and isinstance(func.value, ast.Name) and func.value.id == "json":
        return opened(node.args[0]) if len(node.args) == 1 and not node.keywords else None
    if func.attr in ("read", "readlines") and plain:
        return opened(func.value)
    if func.attr in ("read_text", "read_bytes") and not node.args and literal_keywords(node):
        return opened(func.value, "Path")
    return None


def guarded(node: ast.AST) -> bool:
    """True when the node imports, rebinds, patches or reflects in a way a read cannot allow."""
    if isinstance(node, ast.Import):
        return any(a.name not in ("json", "pathlib") or a.asname for a in node.names)
    if isinstance(node, ast.ImportFrom):
        names = [(a.name, a.asname) for a in node.names]
        return node.module != "pathlib" or node.level != 0 or names != [("Path", None)]
    if isinstance(node, (ast.Global, ast.Nonlocal)):
        return bool(READ_NAMES & set(node.names))
    if isinstance(node, (ast.Attribute, ast.Subscript)) and not isinstance(node.ctx, ast.Load):
        return True
    if isinstance(node, ast.Name):
        stored = not isinstance(node.ctx, ast.Load) and node.id in READ_NAMES
        return stored or node.id in REFLECTION or node.id.startswith("_")
    if isinstance(node, ast.Attribute):
        return node.attr.startswith("_")
    bound = [getattr(node, "name", None), getattr(node, "rest", None), getattr(node, "arg", None)]
    return not isinstance(node, ast.alias) and bool(READ_NAMES & set(bound))


def read_literals(tree: ast.AST) -> set[int]:
    """The ids of the literals that flow only into a read, or none when any node is guarded."""
    nodes = list(ast.walk(tree))
    if any(guarded(node) for node in nodes):
        return set()
    return {id(lit) for node in nodes if (lit := read_of(node)) is not None}
