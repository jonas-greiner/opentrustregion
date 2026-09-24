#!/usr/bin/env python3
# Copyright (C) 2025- Jonas Greiner
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

"""
Check, and optionally fix, the layout of Fortran statements.

Every statement is checked for the column it starts at, which follows from the nesting
of the blocks (modules, procedures, interfaces, derived types, do, if, select, ...) it
sits in, and for how it is broken into continuation lines and spaced. The issues
reported are:

- indent: the statement does not start at the column its block nesting gives
- align: a continuation line does not start where its anchor puts it
- overlong: a line is longer than the column limit
- spacing: operators are not spaced as rule 7 says
- layout: the statement is not broken as the rules below say

The layout of a statement is the best one under these rules:

1. No line is longer than the column limit.
2. The statement uses as few lines as possible, subject only to the first part of rule
   4, and among the layouts that do, breaks each line as late as possible, subject to
   the rest of rule 4.
3. A continuation line is indented relative to the innermost open anchor. Anchors are
   brackets, the assignment operators ``=`` and ``=>``, the ``::`` of a declaration,
   the ``only:`` of a use statement and the ``=`` or ``=>`` of a declaration
   initializer. Breaking directly after an anchor gives a hanging indent of four
   columns, relative to the line a bracket is on or to the start of the statement for
   the other anchors; breaking later aligns with the first thing after the anchor. A
   bracket whose first argument is a bracket group hanging on the same line hangs with
   it. The hanging indent of an initializer in a declaration of several entities is
   instead four columns relative to the names of the entities. Without an open anchor,
   continuation lines are indented by four columns.
4. Bracket groups are kept whole in three tiers. First, before the line count: an array
   constructor, an array section, or the shape of an array in a declaration or an
   allocate statement is not split if it fits whole on a line, and neither are the
   contents of a group holding a single argument or index hung after its opening
   bracket. Second, among the layouts with the fewest lines: a group containing no
   other bracket group apart from grouping parentheses is not split if it fits whole,
   and no group that fits whole has its contents hung after its opening bracket. A
   group fits whole if it fits where it stands, where a break directly before it would
   put it, or on a continuation line with a hanging indent. Where no layout keeps a
   tier, the one breaking it the fewest times is used. Third, the layouts splitting the
   fewest operands are preferred, where an operand is split by a break after an
   operator or comma binding tighter than another one in the same bracket group, and
   counts once for every such group. So a line is broken at ``.or.`` rather than inside
   the comparison next to it, and at a comma rather than inside an argument. An
   operator or comma also binds tighter than any operator or comma outside its bracket
   group, and so does a break directly after an opening bracket, so that a call is kept
   whole where breaking outside it costs no line. The condition of a single-line ``if``
   is looser than the statement it controls, and the commas between the entities of a
   declaration are looser than an initializer.
5. A trailing ``then`` that does not fit takes a hanging indent of four columns like
   any other continuation line. Breaking before it is looser than any break in the
   condition.
6. An array constructor given to ``reshape`` that is written over several lines is
   written one matrix row per line when every row fits. The row length is the first
   extent of a literal shape, or else the number of items per line the constructor is
   already written with. When the rows do not fit, the constructor follows the normal
   rules.
7. Binary ``+``, ``-``, ``*`` and ``/`` have a blank on either side, ``**`` has none,
   and ``//`` is left as written.
8. String literals joined by ``//`` are merged into one, and a string literal is only
   split after a blank between two words, into ``"... "// &`` and ``"..."`` on the next
   line. A literal that fits whole on a continuation line is kept whole like an
   innermost bracket group under rule 4; one that has to be split anyway fills lines
   like any other text, so it may start on the line before.

Statements with comments, ``;`` or a continuation line starting with ``&`` are only
checked for continuation indentation, line length and operator spacing, and are never
reflowed; only their indentation is fixed. Comment lines are not checked, but move with
the statement they are in when its indentation is fixed.

By default only statements touched by the uncommitted changes relative to a git
reference are checked, so that the check can be run after every edit without reporting
the whole code base. Files that git does not track count as changed throughout. Use
--all to check whole files. Run from a subdirectory, only the files below it are
checked.
"""

import argparse
import re
import subprocess
import sys
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

LIMIT = 88
HANG = 4


# lines and statements


def code_part(line: str) -> str:
    """
    this function returns the line with string literals blanked and a trailing comment
    removed, so that column positions stay valid for the original line
    """
    out = []
    quote = None
    for ch in line:
        if quote:
            out.append(ch if ch == quote else " ")
            if ch == quote:
                quote = None
        elif ch in "'\"":
            quote = ch
            out.append(ch)
        elif ch == "!":
            break
        else:
            out.append(ch)
    return "".join(out).rstrip()


def is_code(line: str) -> bool:
    """
    this function returns whether a line carries code, as opposed to being blank, a
    comment or a preprocessor directive
    """
    stripped = line.lstrip()
    return bool(stripped) and not stripped.startswith(("!", "#"))


def indent_of(line: str) -> int:
    """
    this function returns the column the text of a line starts at
    """
    return len(line) - len(line.lstrip())


@dataclass
class Statement:
    """
    this class holds the lines a statement spans, together with those of them that
    carry code
    """

    first: int  # index of the first line
    last: int  # index of the last line, inclusive
    code: List[int] = field(default_factory=list)  # indices of the code lines


def statements(lines: List[str]) -> List[Statement]:
    """
    this function groups code lines into statements, where comment and blank lines
    between the lines of a continued statement belong to it
    """
    result = []
    current = None
    for i, line in enumerate(lines):
        if not is_code(line):
            continue
        if current is None:
            current = Statement(i, i)
        current.code.append(i)
        current.last = i
        if not code_part(line).endswith("&"):
            result.append(current)
            current = None
    if current is not None:
        result.append(current)
    return result


@dataclass
class Joined:
    """
    this class holds a statement joined onto one line, together with where it is
    currently broken and whether it may be reflowed
    """

    text: str  # the statement on one line
    code: str  # the same with string literals blanked
    breaks: List[int]  # positions the statement is currently broken at
    reflowable: bool


def join_statement(lines: List[str], stmt: Statement) -> Joined:
    """
    this function joins the lines of a statement into one, recording where it is
    currently broken
    """
    text, code, breaks = "", "", []
    reflowable = stmt.last - stmt.first + 1 == len(stmt.code)
    for k, i in enumerate(stmt.code):
        line = lines[i]
        blanked = code_part(line)
        if line[len(blanked) :].strip():
            reflowable = False
        start = indent_of(line)
        seg_code, seg_text = blanked[start:], line[start : len(blanked)]
        if k < len(stmt.code) - 1:
            seg_code = seg_code[:-1].rstrip()
            seg_text = seg_text[: len(seg_code)]
        if seg_code.startswith("&"):
            reflowable = False
        if k > 0:
            # no blank after an opening bracket, before a closing one, or after a
            # concatenation written without blanks
            tight = code.endswith(("(", "[")) or seg_code.startswith((")", "]"))
            tight |= code.endswith("//") and not code.endswith(" //")
            text, code = text + ("" if tight else " "), code + ("" if tight else " ")
        text, code = text + seg_text, code + seg_code
        if k < len(stmt.code) - 1:
            breaks.append(len(code))
    if ";" in code:
        reflowable = False
    return Joined(text, code, breaks, reflowable)


# block indentation


PROCEDURE = re.compile(
    r"((pure|impure|elemental|recursive|non_recursive|module)\s+)*"
    r"((integer|real|logical|complex|character|type|class)\s*(\([^)]*\))?\s+)?"
    r"((pure|impure|elemental|recursive|non_recursive|module)\s+)*"
    r"(subroutine|function)\s+\w+"
)


def block_role(code: str) -> Optional[str]:
    """
    this function returns whether a statement closes a block, continues one (else,
    case, contains, ...), opens one, or none of these
    """
    c = code.strip().lower()
    c = re.sub(r"^\w+\s*:\s*(?=(do|if|select|block|associate)\b)", "", c)
    if re.match(
        r"end(\s|$)|end(do|if|select|where|forall|associate|block|interface|type|"
        r"module|subroutine|function|program)\b",
        c,
    ):
        return "close"
    if re.match(
        r"(else|elseif|elsewhere|case|contains)\b|class\s+(is|default)\b|type\s+is\b",
        c,
    ):
        return "middle"
    if re.match(r"module\s+procedure\b", c):
        return None
    if re.match(r"(module|submodule|program)\b", c):
        return "open"
    if re.match(r"(abstract\s+)?interface\b", c) or PROCEDURE.match(c):
        return "open"
    if re.match(r"type\s*(,[^:]*)?::", c) or re.match(r"type\s+\w+\s*$", c):
        return "open"
    if re.match(r"do(\s|$)|block\s*$|critical\s*$|associate\s*\(|enum\s*,", c):
        return "open"
    if re.match(r"select\s+(case|type|rank)\b", c):
        return "open"
    if re.match(r"if\s*\(", c) and re.search(r"\bthen$", c):
        return "open"
    if re.match(r"(where|forall)\s*\(", c):
        depth = 0
        for col, ch in enumerate(c):
            if ch == "(":
                depth += 1
            elif ch == ")":
                depth -= 1
                if depth == 0:
                    return "open" if not c[col + 1 :].strip() else None
    return None


def expected_indents(lines: List[str], stmts: List[Statement]) -> List[int]:
    """
    this function returns the column every statement is expected to start at
    """
    level, result = 0, []
    for stmt in stmts:
        role = block_role(join_statement(lines, stmt).code)
        if role == "close":
            level = max(level - 1, 0)
            result.append(HANG * level)
        elif role == "middle":
            result.append(HANG * max(level - 1, 0))
        else:
            result.append(HANG * level)
            if role == "open":
                level += 1
    return result


# rewrites of a joined statement


def literals(text: str) -> List[Tuple[int, int]]:
    """
    this function returns the first and last position of every string literal in a
    statement
    """
    spans, i = [], 0
    while i < len(text):
        if text[i] in "'\"":
            quote, j = text[i], i + 1
            while j < len(text):
                if text[j] == quote:
                    if text[j + 1 : j + 2] != quote:
                        break
                    j += 1
                j += 1
            if j < len(text):
                spans.append((i, j))
            i = j
        i += 1
    return spans


def merge_literals(joined: Joined) -> Joined:
    """
    this function returns the statement with every concatenation of two string literals
    written as one literal
    """
    text, breaks = joined.text, list(joined.breaks)
    merged = True
    while merged:
        merged = False
        spans = literals(text)
        for (_, a_end), (b_start, _) in zip(spans, spans[1:]):
            same = text[a_end] == text[b_start]
            if same and re.fullmatch(r"\s*//\s*", text[a_end + 1 : b_start]):
                cut = b_start + 1 - a_end
                breaks = [
                    b if b <= a_end else b - cut
                    for b in breaks
                    if not a_end < b <= b_start
                ]
                text = text[:a_end] + text[b_start + 1 :]
                merged = True
                break
    code = code_part(text)
    if len(code) != len(text):
        return joined
    return Joined(text, code, breaks, joined.reflowable)


def next_content(code: str, i: int) -> int:
    """
    this function returns the first position from i on that is not a blank
    """
    while i < len(code) and code[i] == " ":
        i += 1
    return i


def operand_before(code: str, i: int) -> bool:
    """
    this function returns whether the code before position i, ignoring blanks, ends an
    operand, which makes an operator at i binary
    """
    j = i - 1
    while j >= 0 and code[j] == " ":
        j -= 1
    if j < 0:
        return False
    ch = code[j]
    # the sign of an exponent such as 1e-3 is not an operator
    if ch in "eEdD" and j > 0 and code[j - 1] in "0123456789.":
        k = j - 1
        while k >= 0 and code[k] in "0123456789.":
            k -= 1
        if k < 0 or not (code[k].isalnum() or code[k] == "_"):
            return False
    # a dotted operator such as .gt. ends in a period but is no operand
    if ch == "." and re.search(
        r"\.(and|or|eqv|neqv|not|eq|ne|lt|le|gt|ge)\.$", code[: j + 1], re.IGNORECASE
    ):
        return False
    return ch.isalnum() or ch in "_)]\"'."


def space_operators(joined: Joined) -> Joined:
    """
    this function returns the statement with a blank on either side of every binary +,
    -, * and /, and none around **
    """
    code = joined.code
    if re.match(r"\s*(real|integer|logical|complex|character)\s*\*", code, re.I):
        return joined
    inserts, removals = [], []
    for i, ch in enumerate(code):
        if ch not in "+-*/":
            continue
        if code.startswith("**", i) and operand_before(code, i):
            j = i - 1
            while j >= 0 and code[j] == " ":
                removals.append(j)
                j -= 1
            j = i + 2
            while j < len(code) and code[j] == " ":
                removals.append(j)
                j += 1
            continue
        if ch == "*" and "*" in (code[i - 1 : i], code[i + 1 : i + 2]):
            continue
        # the default format of print and read
        if ch == "*" and re.search(r"\b(print|read)\s*$", code[:i], re.I):
            continue
        if ch == "/" and (
            code[i + 1 : i + 2] in ("/", "=", ")") or code[i - 1 : i] in ("/", "(")
        ):
            continue
        if not operand_before(code, i):
            continue
        if code[i - 1] != " ":
            inserts.append(i)
        if code[i + 1 : i + 2] not in (" ", ""):
            inserts.append(i + 1)
    if not inserts and not removals:
        return joined
    text, new_code = joined.text, code
    edits = [(pos, " ") for pos in inserts] + [(pos, "") for pos in removals]
    for pos, blank in sorted(edits, reverse=True):
        end = pos if blank else pos + 1
        text = text[:pos] + blank + text[end:]
        new_code = new_code[:pos] + blank + new_code[end:]
    breaks = [
        b + sum(pos < b for pos in inserts) - sum(pos < b for pos in removals)
        for b in joined.breaks
    ]
    return Joined(text, new_code, breaks, joined.reflowable)


# the structure of a joined statement


DOTTED = re.compile(r"\.(and|or|eqv|neqv|eq|ne|lt|le|gt|ge)\.", re.IGNORECASE)

# how tightly a comma or binary operator binds its operands
PRECEDENCE = {",": 0, ".eqv.": 1, ".neqv.": 1, ".or.": 2, ".and.": 3, "//": 6}
PRECEDENCE.update({op: 5 for op in ("==", "/=", "<", "<=", ">", ">=")})
PRECEDENCE.update({f".{op}.": 5 for op in ("eq", "ne", "lt", "le", "gt", "ge")})
PRECEDENCE.update({"+": 7, "-": 7, "*": 8, "/": 8, "**": 9})


@dataclass
class Structure:
    """
    this class holds the anchors, break candidates and bracket groups of a joined
    statement, which its layouts are chosen from
    """

    # position -> events at that position: ("open", content, start, close, lead, leaf,
    # array, single) for an opening bracket, ("close",) for a closing one, ("anchor",
    # content, kind) for any other anchor and ("entity_end",) for a comma ending an
    # entity of a declaration
    events: Dict[int, List[tuple]]
    candidates: List[int]  # positions a line may be broken at, ascending
    inner: Dict[int, int]  # breaks -> number of groups with a looser break in them
    strings: Dict[int, str]  # breaks inside a string literal -> its quote character
    lengths: Dict[int, int]  # breaks inside a string literal -> length of the literal
    multi_entity: bool  # whether this is a declaration of several entities
    matrices: List[Tuple[int, int, List[int]]]  # (open, close, row ends)


def structure(code: str, text: str, breaks: List[int]) -> Structure:
    """
    this function finds the anchors and break candidates of a joined statement
    """
    events: Dict[int, List[tuple]] = {}
    candidates: Set[int] = set()

    def add(i: int, event: tuple) -> None:
        events.setdefault(i, []).append(event)

    lower = code.lower()
    end = len(code.rstrip())
    is_use = re.match(r"\s*use\b", lower) is not None

    # matching brackets
    match, stack = {}, []
    for i, ch in enumerate(code):
        if ch in "([":
            stack.append(i)
        elif ch in ")]" and stack:
            match[stack.pop()] = i

    # a declaration has a :: outside brackets
    depth, decl_at = 0, None
    for i, ch in enumerate(code):
        if ch in "([":
            depth += 1
        elif ch in ")]":
            depth -= 1
        elif depth == 0 and code.startswith("::", i) and not is_use:
            decl_at = i
            break

    # parentheses that only group an expression, as opposed to those of a call, an
    # array reference or a statement keyword such as if
    def grouping(o: int) -> bool:
        j = o - 1
        while j >= 0 and code[j] == " ":
            j -= 1
        return code[o] == "(" and (j < 0 or not (code[j].isalnum() or code[j] in "_%"))

    # how tightly the operator or comma every break follows binds, in the bracket group
    # it is in and in the groups around it, where it binds tighter than any operator or
    # comma outside its bracket group
    binds: Dict[int, List[Tuple[int, int]]] = {}
    opening: Dict[int, List[Tuple[int, int]]] = {}
    groups = [-1]

    def bind(rank: int) -> List[Tuple[int, int]]:
        entries = [(rank, groups[-1])]
        for k in range(len(groups) - 1, 0, -1):
            rank += 10
            entries.append((rank, groups[k - 1]))
        return entries

    depth, assigned, in_list, entities = 0, None, False, 0
    initializers = []
    i = 0
    while i < end:
        ch = code[i]
        if ch in "([":
            groups.append(i)
            # the group a bracket belongs to starts at the name before it
            start = i
            while start > 0 and (code[start - 1].isalnum() or code[start - 1] in "_%"):
                start -= 1
            add(i, ("open", next_content(code, i + 1), start, match.get(i, end)))
            depth += 1
            candidates.add(i + 1)
            # a break directly after the bracket is inside its group, and so binds
            # tighter than anything in the groups around it
            opening[i + 1] = bind(PRECEDENCE[","])[1:]
            i += 1
            continue
        if ch in ")]":
            add(i, ("close",))
            depth -= 1
            if len(groups) > 1:
                groups.pop()
            i += 1
            continue
        if ch == ",":
            candidates.add(i + 1)
            binds[i + 1] = bind(PRECEDENCE[","])
            if in_list and depth == 0:
                add(i, ("entity_end",))
                entities += 1
            i += 1
            continue
        if decl_at is not None and i == decl_at:
            add(i + 1, ("anchor", next_content(code, i + 2), "list"))
            candidates.add(i + 2)
            in_list, entities = True, 1
            i += 2
            continue
        if is_use and depth == 0 and lower.startswith("only", i):
            colon = next_content(code, i + 4)
            if code[colon : colon + 1] == ":":
                add(colon, ("anchor", next_content(code, colon + 1), "list"))
                candidates.add(colon + 1)
                i = colon + 1
                continue
        if ch == "=" and depth == 0 and not is_use:
            prev = code[i - 1] if i > 0 else " "
            nxt = code[i + 1] if i + 1 < len(code) else " "
            if nxt != "=" and prev not in "=/<>":
                op_end = i + 2 if nxt == ">" else i + 1
                if in_list:
                    add(op_end - 1, ("anchor", next_content(code, op_end), "init"))
                    candidates.add(op_end)
                    initializers.append(op_end)
                elif not assigned and decl_at is None:
                    add(op_end - 1, ("anchor", next_content(code, op_end), "assign"))
                    candidates.add(op_end)
                    assigned = op_end
                i = op_end
                continue
        # binary operators
        op = None
        m = DOTTED.match(code, i)
        if m:
            op = m.end()
        elif code.startswith(("//", "**", "==", "/=", "<=", ">="), i):
            if operand_before(code, i):
                op = i + 2
        elif ch in "+-*/<>" and operand_before(code, i):
            if not (ch == "/" and code[i + 1 : i + 2] == ")"):
                op = i + 1
        if op is not None:
            candidates.add(op)
            binds[op] = bind(PRECEDENCE[lower[i:op].strip()])
            i = op
            continue
        i += 1

    # in a declaration of several entities, the initializer of an entity binds tighter
    # than the commas separating the entities
    if entities > 1:
        for b in initializers:
            binds[b] = [(1, -1)]

    # a break after the condition of a single-line if or where statement, which is
    # looser than any break in the statement it controls
    m = re.match(r"\s*(if|where)\s*\(", lower)
    if m and match.get(m.end() - 1) is not None:
        after = match[m.end() - 1] + 1
        if not re.match(r"\s*then\s*$", lower[after:]):
            candidates.add(after)
            binds[after] = [(-2, -1)]
            # the assignment it controls, which binds tighter than the condition
            if assigned is not None and assigned > after:
                binds[assigned] = [(-1, -1)]

    # breaks inside string literals, after a blank between two words
    strings: Dict[int, str] = {}
    lengths: Dict[int, int] = {}
    for first, last in literals(text):
        for b in range(first + 3, last):
            if text[b - 1] == " " and text[b] != " " and text[b - 2] != " ":
                strings[b] = text[first]
                lengths[b] = last - first + 1

    # a break before a lone then, which is looser than any break in the condition, and
    # before the result or bind of a procedure
    m = re.search(r"\)\s*(then)\s*$", lower)
    if m:
        candidates.add(m.start(1))
        binds[m.start(1)] = [(-2, -1)]
    for kw in ("result", "bind"):
        for m in re.finditer(rf"\)\s*({kw})\s*\(", lower):
            candidates.add(m.start(1))

    # a break that splits an operand of a looser operator or comma in the same group
    loosest: Dict[int, int] = {}
    for entries in binds.values():
        for rank, group in entries:
            loosest[group] = min(rank, loosest.get(group, rank))
    inner = {
        b: sum(rank > loosest.get(group, rank) for rank, group in entries)
        for b, entries in list(binds.items()) + list(opening.items())
    }

    candidates = {b for b in candidates if 0 < b < end and next_content(code, b) < end}
    candidates.update(strings)

    # an array constructor, an array section, or the shape of an array in a declaration
    # or an allocate statement
    allocating = re.match(r"\s*(de)?allocate\s*\(", lower)

    def array_like(o: int, c: int, enclosing: List[int]) -> bool:
        if code[o] == "[":
            return True
        depth = 0
        for j in range(o + 1, c):
            if code[j] in "([":
                depth += 1
            elif code[j] in ")]":
                depth -= 1
            elif depth == 0 and code[j] == ":" and ":" not in code[j - 1 : j + 2 : 2]:
                return True
        if decl_at is not None and o > decl_at and not enclosing:
            return True
        return bool(allocating) and len(enclosing) == 1

    # a bracket group holding a single argument, index or expression
    def single(o: int, c: int) -> bool:
        depth = 0
        for j in range(o + 1, c):
            if code[j] in "([":
                depth += 1
            elif code[j] in ")]":
                depth -= 1
            elif depth == 0 and code[j] == ",":
                return False
        return True

    # a bracket group can be moved to the next line whole, together with what joins it
    # to the last possible break before it (such as .not.), when that break lies inside
    # the bracket group enclosing it
    ordered = sorted(candidates - set(strings))
    for i, evs in events.items():
        for k, ev in enumerate(evs):
            if ev[0] == "open":
                before = [b for b in ordered if next_content(code, b) <= ev[2]]
                lead = next_content(code, before[-1]) if before else None
                enclosing = [o for o, c in match.items() if o < ev[2] and c > ev[3]]
                if lead is not None and enclosing and lead <= max(enclosing):
                    lead = None
                leaf = not any(i < o < ev[3] and not grouping(o) for o in match)
                array = array_like(i, ev[3], enclosing)
                evs[k] = ev + (lead, leaf, array, single(i, ev[3]))

    # array constructors passed to reshape, laid out one matrix row per line
    matrices = []
    for m in re.finditer(r"\breshape\s*\(\s*\[", lower):
        open_at = m.end() - 1
        close_at = match.get(open_at)
        if close_at is None:
            continue
        depth, commas = 0, []
        for i in range(open_at + 1, close_at):
            if code[i] in "([":
                depth += 1
            elif code[i] in ")]":
                depth -= 1
            elif code[i] == "," and depth == 0:
                commas.append(i)
        n_items = len(commas) + 1
        row = None
        written = [b for b in breaks if open_at < b < close_at]
        shape = re.match(r"\s*,\s*\[\s*(\d+)\s*,", code[close_at + 1 :])
        if not written:
            pass
        elif shape:
            row = int(shape.group(1))
        else:
            ends = [next_content(code, c + 1) for c in commas]
            if all(next_content(code, b) in ends for b in written):
                bounds = [ends.index(next_content(code, b)) + 1 for b in written]
                bounds.append(n_items)
                sizes = [b - a for a, b in zip([0] + bounds[:-1], bounds)]
                if len(set(sizes)) == 1:
                    row = sizes[0]
        if row and 1 < row < n_items and n_items % row == 0:
            row_ends = [commas[k * row - 1] + 1 for k in range(1, n_items // row)]
            matrices.append((open_at, close_at, row_ends))

    return Structure(
        events, sorted(candidates), inner, strings, lengths, entities > 1, matrices
    )


# the layout search


# flags of an open bracket group under rule 4: it fits whole on a line, contains no
# other bracket group apart from grouping parentheses, is array-like, or holds a single
# argument or index
FITS, LEAF, ARRAY, SINGLE = 1, 2, 4, 8

# an open anchor: its kind, the column of the first thing after it (None when that
# hangs), its hanging indentation, its rule 4 flags and the position of the first thing
# after it
Anchor = Tuple[str, Optional[int], int, int, int]

# the cost of a layout in the order it is minimized (array-like groups split, lines,
# other groups split, operands split), and its breaks, negated so that the latest
# breaks sort first
Score = Tuple[int, int, int, int, Tuple[int, ...]]


class Layout:
    """
    this class holds the layouts of one joined statement starting at a given column
    """

    def __init__(self, joined: Joined, base: int, matrix: bool):
        self.j = joined
        self.base = base
        self.s = structure(joined.code, joined.text, joined.breaks)
        self.end = len(joined.code.rstrip())
        self.forced: List[int] = []
        self.inside: List[Tuple[int, int]] = []
        if matrix:
            for open_at, close_at, row_ends in self.s.matrices:
                self.forced += row_ends
                self.inside.append((open_at, close_at))
        self.forced.sort()
        # the start of the bracket group, name included, whose contents begin at a
        # position
        self.starts = {
            ev[1]: ev[2]
            for evs in self.s.events.values()
            for ev in evs
            if ev[0] == "open"
        }
        # a string literal that fits whole on a continuation line is split only when
        # that saves a line
        self.whole_strings = {
            b
            for b, length in self.s.lengths.items()
            if base + HANG + length + len(", &") <= LIMIT
        }
        self.candidates = [
            b
            for b in self.s.candidates
            if b in self.forced or not any(a < b < c for a, c in self.inside)
        ]

    def scan(
        self, s: int, b: int, indent: int, stack: Tuple[Anchor, ...]
    ) -> Tuple[Anchor, ...]:
        """
        this function returns the open anchors after laying out positions s to b on a
        line starting at the given indentation
        """

        # a line continuing a string literal starts with its quote
        pre = int(s in self.s.strings)

        def col(i: int) -> int:
            return indent + pre + i - s

        anchors = list(stack)
        for i in range(s, b):
            for ev in self.s.events.get(i, ()):
                if ev[0] == "open":
                    _, content, start, close, lead, leaf, array, one = ev
                    fits = False
                    # a group holding matrix rows is split by them
                    rows = any(i <= a and c <= close for a, c in self.inside)
                    if lead is not None and lead >= s and not rows:
                        # where the group would start if moved to the next line
                        if lead == s:
                            at = indent
                        elif anchors:
                            _, anchor, hang, _, top_content = anchors[-1]
                            if anchor is None or lead <= top_content:
                                at = hang
                            else:
                                at = anchor
                        else:
                            at = self.base + HANG
                        # whether it fits there or on a continuation line with a
                        # hanging indent, which an earlier break can give it
                        at = min(at, self.base + HANG)
                        tail = 2 if close + 1 < self.end else 0
                        fits = at + close - lead + 1 + tail <= LIMIT
                    # an innermost group that fits whole where it stands is not split
                    # either
                    tail = 2 if close + 1 < self.end else 0
                    here = col(start) + close - start + 1 + tail <= LIMIT
                    fits = fits or (leaf and not rows and here)
                    whole = 0
                    if fits:
                        whole = FITS | LEAF * leaf | ARRAY * array | SINGLE * one
                    anchor = col(content) if content < b else None
                    if anchor is None:
                        # a bracket whose first argument hangs hangs with it
                        pos = start
                        for k in range(len(anchors) - 1, -1, -1):
                            kind, outer, hang, outer_whole, outer_content = anchors[k]
                            if kind != "bracket" or outer_content != pos:
                                break
                            anchors[k] = (kind, None, hang, outer_whole, outer_content)
                            pos = self.starts[outer_content]
                    anchors.append(("bracket", anchor, indent + HANG, whole, content))
                elif ev[0] == "close":
                    if anchors and anchors[-1][0] == "bracket":
                        anchors.pop()
                elif ev[0] == "anchor":
                    _, content, kind = ev
                    names = next((a for a in anchors if a[0] == "list"), None)
                    if kind == "init" and self.s.multi_entity and names is not None:
                        hang = (names[2] if names[1] is None else names[1]) + HANG
                    else:
                        hang = self.base + HANG
                    anchor = col(content) if content < b else None
                    anchors.append((kind, anchor, hang, 0, content))
                elif ev[0] == "entity_end":
                    if anchors and anchors[-1][0] == "init":
                        anchors.pop()
        return tuple(anchors)

    def next_indent(self, b: int, stack: Tuple[Anchor, ...]) -> int:
        """
        this function returns the indentation of the line following a break
        """
        rest = self.j.code[next_content(self.j.code, b) :].lower()
        if re.match(r"then\s*$|(result|bind)\s*\(", rest):
            return self.base + HANG
        if stack:
            _, anchor, hang, _, _ = stack[-1]
            return hang if anchor is None else anchor
        return self.base + HANG

    def best(self) -> Optional[List[int]]:
        """
        this function returns the break positions of the best layout, or None when no
        layout keeps within the column limit
        """
        text, end = self.j.text, self.end

        @lru_cache(maxsize=None)
        def search(s: int, indent: int, stack: Tuple[Anchor, ...]) -> Optional[Score]:
            options: List[Score] = []
            forced = next((f for f in self.forced if f > s), None)
            start = indent + int(s in self.s.strings)
            if forced is None and start + len(text[s:end].rstrip()) <= LIMIT:
                options.append((0, 1, 0, 0, ()))
            for b in reversed(self.candidates):
                if b <= s or (forced is not None and b > forced):
                    continue
                if b in self.s.strings:
                    # the literal is closed and concatenated with its continuation
                    width = start + len(text[s:b]) + len('"// &')
                else:
                    width = start + len(text[s:b].rstrip()) + len(" &")
                if width > LIMIT:
                    continue
                new_stack = self.scan(s, b, indent, stack)
                # splitting an array-like group that fits whole, and splitting an
                # innermost group or hanging the contents of any group that fits whole
                arrays = groups = 0
                for _, anchor, _, whole, _ in new_stack:
                    hung = anchor is None
                    if whole & ARRAY or (whole & SINGLE and hung):
                        arrays += 1
                    elif whole & LEAF or (whole and hung):
                        groups += 1
                split = self.s.inner.get(b, 0)
                groups += b in self.whole_strings
                n = b if b in self.s.strings else next_content(self.j.code, b)
                sub = search(n, self.next_indent(b, new_stack), new_stack)
                if sub is not None:
                    options.append(
                        (
                            sub[0] + arrays,
                            sub[1] + 1,
                            sub[2] + groups,
                            sub[3] + split,
                            (-b,) + sub[4],
                        )
                    )
            # fewest array-like groups split, then fewest lines, then fewest other
            # groups split against rule 4, then fewest split operands, then the latest
            # first break, the latest second break, ...
            return min(options) if options else None

        result = search(next_content(self.j.code, 0), self.base, ())
        return None if result is None else [-b for b in result[4]]

    def render(self, breaks: List[int]) -> List[str]:
        """
        this function returns the lines of the statement broken at the given positions
        """
        lines, indent = [], self.base
        stack: Tuple[Anchor, ...] = ()
        s = next_content(self.j.code, 0)
        for b in breaks:
            line = " " * indent + self.s.strings.get(s, "")
            if b in self.s.strings:
                line += self.j.text[s:b] + self.s.strings[b] + "// &"
            else:
                line += self.j.text[s:b].rstrip() + " &"
            lines.append(line)
            stack = self.scan(s, b, indent, stack)
            indent = self.next_indent(b, stack)
            s = b if b in self.s.strings else next_content(self.j.code, b)
        start = " " * indent + self.s.strings.get(s, "")
        lines.append(start + self.j.text[s : self.end].rstrip())
        return lines


def best_layout(joined: Joined, base: int) -> Optional[List[str]]:
    """
    this function returns the best layout of a statement, dropping the matrix rows when
    no layout keeps them, and keeping string literals as they are written when merging
    them leaves no layout at all
    """
    for statement in (merge_literals(joined), joined):
        statement = space_operators(statement)
        for matrix in (True, False):
            layout = Layout(statement, base, matrix)
            if matrix and not layout.forced:
                continue
            breaks = layout.best()
            if breaks is not None:
                return layout.render(breaks)
    return None


# checking and fixing


@dataclass
class Issue:
    """
    this class holds an issue of a statement, together with what fixes it
    """

    line: int
    kind: str
    message: str
    stmt: Statement
    replacement: Optional[List[str]] = None  # new lines of the whole statement
    indent: Optional[int] = None  # the column the line should start at


def in_scope(stmt: Statement, touched: Optional[Set[int]]) -> bool:
    """
    this function returns whether a statement is to be checked, which it is when it
    contains a touched line or when every line counts as touched
    """
    return touched is None or bool(
        touched.intersection(range(stmt.first, stmt.last + 1))
    )


def analyse(lines: List[str], touched: Optional[Set[int]]) -> List[Issue]:
    """
    this function returns the indentation and layout issues of a file, restricted to
    statements containing a touched line when a set of touched lines is given
    """
    issues = []
    stmts = statements(lines)
    for stmt, base in zip(stmts, expected_indents(lines, stmts)):
        if not in_scope(stmt, touched):
            continue
        have = indent_of(lines[stmt.first])
        if have != base:
            message = f"starts at {have}, expected {base}"
            issues.append(Issue(stmt.first, "indent", message, stmt, indent=base))
            continue
        joined = join_statement(lines, stmt)
        if not joined.reflowable:
            # keep the breaks, only check where the lines start and end
            layout = Layout(joined, base, False)
            current = [lines[i].rstrip() for i in stmt.code]
            want = layout.render(joined.breaks)
            for i, have_line, want_line in zip(stmt.code, current, want):
                if indent_of(have_line) != indent_of(want_line):
                    column = indent_of(want_line)
                    message = f"indent {indent_of(have_line)}, expected {column}"
                    issues.append(Issue(i, "align", message, stmt, indent=column))
                elif len(have_line) > LIMIT:
                    message = f"{len(have_line)} columns, statement not reflowable"
                    issues.append(Issue(i, "overlong", message, stmt))
            if space_operators(joined).code != joined.code:
                message = "operator spacing, statement not reflowable"
                issues.append(Issue(stmt.first, "spacing", message, stmt))
            continue
        current = [lines[i].rstrip() for i in range(stmt.first, stmt.last + 1)]
        best = best_layout(joined, base)
        if best is None:
            for i in stmt.code:
                if len(lines[i].rstrip()) > LIMIT:
                    message = f"{len(lines[i].rstrip())} columns, no layout fits"
                    issues.append(Issue(i, "overlong", message, stmt))
            continue
        if best != current:
            if any(len(line) > LIMIT for line in current):
                kind = "overlong"
            elif [line.strip() for line in best] == [line.strip() for line in current]:
                kind = "align"
            elif [re.sub(r"\s", "", line) for line in best] == [
                re.sub(r"\s", "", line) for line in current
            ]:
                kind = "spacing"
            else:
                kind = "layout"
            message = f"{len(current)} line(s), best layout has {len(best)}"
            issues.append(Issue(stmt.first, kind, message, stmt, best))
    return issues


def fix(lines: List[str], touched: Optional[Set[int]]) -> Tuple[List[str], int]:
    """
    this function returns the lines with every fixable issue fixed, and the number of
    fixes; the indentation of statements is fixed first, and alone, since the layout of
    a statement depends on the column it starts at
    """
    changed = 0
    for issue in analyse(lines, touched):
        if issue.kind != "indent" or issue.indent is None:
            continue
        stmt = issue.stmt
        shift = issue.indent - indent_of(lines[stmt.first])
        for i in range(stmt.first, stmt.last + 1):
            if lines[i].strip():
                line = lines[i]
                lines[i] = " " * max(indent_of(line) + shift, 0) + line.lstrip()
        changed += 1
    if changed:
        return lines, changed
    for issue in reversed(analyse(lines, touched)):
        stmt = issue.stmt
        if issue.replacement is not None:
            lines = lines[: stmt.first] + issue.replacement + lines[stmt.last + 1 :]
            changed += 1
        elif issue.kind == "align" and issue.indent is not None:
            lines[issue.line] = " " * issue.indent + lines[issue.line].lstrip()
            changed += 1
    return lines, changed


# command line


def git(*args: str) -> str:
    """
    this function returns the output of a git command, stopping with its error when it
    fails, so that a failing command is never mistaken for one reporting no changes
    """
    result = subprocess.run(
        ["git", "--no-pager", *args], capture_output=True, text=True
    )
    if result.returncode != 0:
        # the first line holds the reason, what follows is usage help
        reason = (result.stderr.strip().splitlines() or ["unknown error"])[0]
        sys.exit(f"git {' '.join(args)} failed: {reason}")
    return result.stdout


def touched_lines(path: Path, ref: str) -> Optional[Set[int]]:
    """
    this function returns the 0-based indices of lines added or changed relative to a
    git reference, or None when the file is not tracked, in which case every line
    counts as touched
    """
    tracked = subprocess.run(
        ["git", "ls-files", "--error-unmatch", str(path)], capture_output=True
    )
    if tracked.returncode != 0:
        return None
    diff = git("diff", "--no-color", "--no-ext-diff", "-U0", ref, "--", str(path))
    result: Set[int] = set()
    for hunk in re.finditer(r"^@@ -\S+ \+(\d+)(?:,(\d+))? @@", diff, re.MULTILINE):
        start, count = int(hunk.group(1)), int(hunk.group(2) or 1)
        if count == 0:
            # a hunk that only deletes lines touches the lines on either side of it
            result.update(i for i in (start - 1, start) if i >= 0)
        else:
            result.update(range(start - 1, start - 1 + count))
    return result


def main() -> int:
    """
    this function checks the files named or found through git, fixes what is reported
    when asked to, and returns 1 when issues remain and 0 otherwise
    """
    parser = argparse.ArgumentParser(description=__doc__.strip().split("\n")[0])
    parser.add_argument("files", nargs="*", type=Path, help="Fortran files to check")
    parser.add_argument("--ref", default="HEAD", help="git reference (default HEAD)")
    parser.add_argument(
        "--all", action="store_true", help="check whole files, not only changes"
    )
    parser.add_argument("--fix", action="store_true", help="fix what is reported")
    parser.add_argument(
        "--show",
        action="store_true",
        help="print the best layout of every statement reported",
    )
    args = parser.parse_args()

    if args.files:
        files = args.files
    elif args.all and args.fix:
        parser.error("--all --fix needs the files to fix to be named explicitly")
    else:
        listings = [["ls-files", "-z"]]
        if not args.all:
            # changed files, and files git does not track yet
            listings = [
                ["diff", "--no-color", "--name-only", "--relative", "-z", args.ref],
                ["ls-files", "-z", "--others", "--exclude-standard"],
            ]
        out = []
        for listing in listings:
            out += git(*listing).split("\0")
        files = [Path(f) for f in out if f.endswith(".f90") and Path(f).exists()]

    remaining = 0
    for path in files:
        lines = path.read_text().split("\n")

        def scope() -> Optional[Set[int]]:
            return None if args.all else touched_lines(path, args.ref)

        if args.fix:
            total, written = 0, "\n".join(lines)
            for _ in range(10):
                lines, changed = fix(lines, scope())
                if not changed:
                    break
                # never overwrite a change saved to the file in the meantime
                if path.read_text() != written:
                    print(f"{path}: changed while fixing, left as it is")
                    lines = path.read_text().split("\n")
                    break
                written = "\n".join(lines)
                path.write_text(written)
                total += changed
            else:
                print(f"{path}: fixes did not settle")
            if total:
                print(f"{path}: {total} fix(es) applied")
        for issue in analyse(lines, scope()):
            print(f"{path}:{issue.line + 1}: {issue.kind}: {issue.message}")
            if args.show and issue.replacement:
                print("\n".join("    | " + line for line in issue.replacement))
            remaining += 1
    return 1 if remaining else 0


if __name__ == "__main__":
    sys.exit(main())
