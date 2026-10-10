#!/usr/bin/env python3
"""Token audit for the current Pi project.

Pi records one JSONL session log per session under
``<agent-dir>/sessions/--<cwd with separators dashed>--/`` and each assistant message
carries a ``usage`` block (input, output, cacheRead, cacheWrite, reasoning). That is the
only trustworthy local source of token accounting for this repository, so this script
reads it instead of guessing.

Usage:
    python3 scripts/token_audit.py                    # table for this working directory
    python3 scripts/token_audit.py --json             # same data, machine readable
    python3 scripts/token_audit.py --save snap.json   # write a snapshot
    python3 scripts/token_audit.py --compare snap.json

It is read-only. It never touches the session logs.

What the numbers mean:
    cacheRead   the whole accumulated context, re-read every turn. It dominates because a
                session's cost is roughly the sum of its per-turn context, i.e. quadratic
                in the turn count. It is cheap per token, but it is the volume to shrink.
    input       context sent without a cache hit — the part that is genuinely new.
    output +    what the model produced.
    reasoning
    tool chars  payload the tools actually returned, estimated at ~4 chars/token. Compare
                it against the totals: when it is ~1% of the total, the conversation is
                re-reading history, not reading files.
"""

from __future__ import annotations

import argparse
import json
import os
import statistics
import sys
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

USAGE_FIELDS = ("input", "output", "cacheRead", "cacheWrite", "reasoning", "totalTokens")
CHARS_PER_TOKEN = 4
SECTIONS = (
    "addendum",
    "skills",
    "project_context",
    "rules",
    "tools",
    "docs",
    "preamble",
    "mcp_servers",
    "cwd",
)


def default_sessions_dir(cwd: Path) -> Path:
    agent_dir = Path(os.environ.get("PI_CODING_AGENT_DIR", Path.home() / ".pi" / "agent"))
    slug = "--" + str(cwd).strip("/").replace("/", "-") + "--"
    return agent_dir / "sessions" / slug


@dataclass
class Session:
    name: str
    turns: int = 0
    usage: Counter = field(default_factory=Counter)
    tool_chars: Counter = field(default_factory=Counter)
    tool_calls: Counter = field(default_factory=Counter)
    contexts: list[int] = field(default_factory=list)
    system_sections: dict[str, int] = field(default_factory=dict)
    size_bytes: int = 0

    @property
    def mean_context(self) -> int:
        return int(statistics.mean(self.contexts)) if self.contexts else 0

    @property
    def peak_context(self) -> int:
        return max(self.contexts) if self.contexts else 0


def load_sessions(sessions_dir: Path, limit: int) -> list[Session]:
    if not sessions_dir.is_dir():
        sys.exit(f"token_audit: no session directory at {sessions_dir}")
    paths = sorted(sessions_dir.glob("*.jsonl"), key=lambda p: p.stat().st_mtime)
    if limit:
        paths = paths[-limit:]
    return [parse_session(p) for p in paths]


def parse_session(path: Path) -> Session:
    session = Session(name=path.name)
    session.size_bytes = path.stat().st_size
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            try:
                record = json.loads(line)
            except ValueError:
                continue
            if record.get("type") != "message":
                continue
            message = record.get("message") or {}
            role = message.get("role")
            if role == "system":
                sections = message.get("sections") or {}
                for name, body in sections.items():
                    size = len(body) if isinstance(body, str) else len(json.dumps(body))
                    # Keep the largest emission: compaction re-emits the prompt several
                    # times per long session and the floor is the number that matters.
                    session.system_sections[name] = max(session.system_sections.get(name, 0), size)
            elif role == "assistant":
                session.turns += 1
                usage = message.get("usage") or {}
                for key in USAGE_FIELDS:
                    session.usage[key] += usage.get(key, 0) or 0
                session.contexts.append(usage.get("totalTokens", 0) or 0)
            elif role == "toolResult":
                name = message.get("toolName", "?")
                session.tool_calls[name] += 1
                session.tool_chars[name] += len(json.dumps(message.get("content")))
    return session


def render(sessions: list[Session], top: int) -> str:
    usage: Counter = Counter()
    tool_chars: Counter = Counter()
    tool_calls: Counter = Counter()
    system_sections: dict[str, int] = {}
    turns = 0
    for session in sessions:
        turns += session.turns
        usage.update(session.usage)
        tool_chars.update(session.tool_chars)
        tool_calls.update(session.tool_calls)
        for name, size in session.system_sections.items():
            system_sections[name] = max(system_sections.get(name, 0), size)

    total = usage["totalTokens"]
    tool_total = sum(tool_chars.values())
    system_total = sum(system_sections.values())
    mean_context = usage["cacheRead"] // turns if turns else 0
    peak_context = max((s.peak_context for s in sessions), default=0)

    lines: list[str] = []
    lines.append(f"sessions: {len(sessions)}   turns: {turns}")
    lines.append(
        f"tokens: total={total:,}  cacheRead={usage['cacheRead']:,} "
        f"input={usage['input']:,}  output={usage['output']:,}  reasoning={usage['reasoning']:,}"
    )
    lines.append(
        f"context: mean/turn={mean_context:,}  peak={peak_context:,}"
        + (
            f"  tool payload≈{tool_total // CHARS_PER_TOKEN:,} tok ({tool_total:,} chars)"
            if tool_total
            else ""
        )
    )
    if total and tool_total:
        lines.append(
            f"  → tool payload is {100 * tool_total / CHARS_PER_TOKEN / total:.2f}% of billed tokens"
        )

    if system_sections:
        rendered = "  ".join(
            f"{name}={system_sections[name] // CHARS_PER_TOKEN:,}"
            for name in SECTIONS
            if name in system_sections
        )
        lines.append(
            f"system prompt floor: ≈{system_total // CHARS_PER_TOKEN:,} tok/turn  [{rendered}]"
        )

    if tool_chars:
        lines.append("top tools by returned payload:")
        for name, chars in tool_chars.most_common(top):
            lines.append(
                f"  {name:20s} n={tool_calls[name]:5d}  chars={chars:12,d}  ≈{chars // CHARS_PER_TOKEN:9,d} tok"
            )

    lines.append("sessions (newest last, ordered by cacheRead):")
    for session in sorted(sessions, key=lambda s: -s.usage["cacheRead"])[:top]:
        lines.append(
            f"  {session.name[:30]:30s} turns={session.turns:4d}  "
            f"cacheRead={session.usage['cacheRead']:13,d}  out={session.usage['output']:9,d}  "
            f"ctxMed={int(statistics.median(session.contexts)) if session.contexts else 0:8,d}"
        )
    return "\n".join(lines)


def as_json(sessions: list[Session]) -> str:
    payload = {
        "sessions": [
            {
                "name": s.name,
                "turns": s.turns,
                "usage": dict(s.usage),
                "toolChars": dict(s.tool_chars),
                "toolCalls": dict(s.tool_calls),
                "meanContext": s.mean_context,
                "peakContext": s.peak_context,
                "systemSections": s.system_sections,
                "sizeBytes": s.size_bytes,
            }
            for s in sessions
        ],
        "totals": {
            "sessions": len(sessions),
            "turns": sum(s.turns for s in sessions),
            "usage": dict(sum((Counter(s.usage) for s in sessions), Counter())),
            "toolChars": dict(sum((Counter(s.tool_chars) for s in sessions), Counter())),
        },
    }
    return json.dumps(payload, indent=2, sort_keys=True)


def compare(previous: dict, current: dict) -> str:
    lines = ["delta vs snapshot:"]
    prev_totals = previous.get("totals", {})
    cur_totals = current["totals"]

    def fmt(delta: int) -> str:
        sign = "+" if delta > 0 else ""
        return f"{sign}{delta:,}"

    prev_usage = prev_totals.get("usage", {})
    cur_usage = cur_totals.get("usage", {})
    for key in ("totalTokens", "cacheRead", "input", "output"):
        before, after = prev_usage.get(key, 0), cur_usage.get(key, 0)
        pct = f" ({100 * (after - before) / before:+.1f}%)" if before else ""
        lines.append(f"  {key:12s} {before:>14,} → {after:>14,}  {fmt(after - before)}{pct}")
    prev_turns, cur_turns = prev_totals.get("turns", 0), cur_totals.get("turns", 0)
    lines.append(
        f"  {'turns':12s} {prev_turns:>14,} → {cur_turns:>14,}  {fmt(cur_turns - prev_turns)}"
    )
    before_per_turn = prev_usage.get("cacheRead", 0) / max(prev_turns, 1)
    after_per_turn = cur_usage.get("cacheRead", 0) / max(cur_turns, 1)
    if before_per_turn:
        delta = 100 * (after_per_turn - before_per_turn) / before_per_turn
        lines.append(
            f"  {'cacheRead/turn':12s} {before_per_turn:>14,.0f} → {after_per_turn:>14,.0f}  ({delta:+.1f}%)"
        )
    lines.append("  (per-turn cacheRead is the metric to watch: it is the re-read floor)")
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--sessions-dir", type=Path, default=None, help="override the Pi session directory"
    )
    parser.add_argument(
        "--cwd", type=Path, default=Path.cwd(), help="project whose sessions to audit"
    )
    parser.add_argument("--top", type=int, default=10, help="rows per ranked section")
    parser.add_argument("--limit", type=int, default=0, help="only the newest N session files")
    parser.add_argument("--json", action="store_true", help="machine-readable output")
    parser.add_argument("--save", type=Path, default=None, help="write a snapshot for --compare")
    parser.add_argument("--compare", type=Path, default=None, help="diff against a saved snapshot")
    args = parser.parse_args()

    sessions_dir = args.sessions_dir or default_sessions_dir(args.cwd)
    sessions = load_sessions(sessions_dir, args.limit)
    payload = json.loads(as_json(sessions))

    if args.compare:
        print(compare(json.loads(args.compare.read_text(encoding="utf-8")), payload))
        return 0
    if args.json:
        print(json.dumps(payload, indent=2, sort_keys=True))
    else:
        print(render(sessions, args.top))
    if args.save:
        args.save.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
        print(f"\nsnapshot written to {args.save}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
