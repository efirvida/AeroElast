#!/usr/bin/env python3
"""Run a Windows executable under the WINE_AppImage inside a pseudo-terminal.

Fixes Intel Fortran runtime errors like
  forrtl: severe (38): error during write, unit -1, file CONOUT$
which occur when the program has no real console (CONOUT$). stdout/stderr are
wired to the pty; stdin is a pipe so console reads can be driven from Python.

Usage: wine-appimage-run-pty.py <windows-program.exe> [args...]

Set WINE_PTY_AUTO_ANSWER (e.g. "n") to answer the Fortran prompt
  Stop execution now? [y]es / [n]o :
automatically; WINE_PTY_PROMPT_RE overrides the prompt regex.
"""
import os
import pty
import re
import select
import signal
import subprocess
import sys
import time

_ESC_SEQ = re.compile(r"\x1b\[[0-9;?]*[A-Za-z]")

AUTO_ANSWER = os.environ.get("WINE_PTY_AUTO_ANSWER")
AUTO_PROMPT = re.compile(os.environ.get(
    "WINE_PTY_PROMPT_RE", r"Stop execution now\? \[y\]es / \[n\]o"))

WINE_DIR = "/scratch/leahk/eduardo.donestevez/AppDir"
LOADER = f"{WINE_DIR}/runtime/compat/lib64/ld-linux-x86-64.so.2"
LIBRARY_PATH = ":".join([
    f"{WINE_DIR}/runtime/compat/lib64",
    f"{WINE_DIR}/runtime/compat/lib/x86_64-linux-gnu",
    f"{WINE_DIR}/lib/x86_64-linux-gnu",
    f"{WINE_DIR}/usr/lib/x86_64-linux-gnu",
    f"{WINE_DIR}/opt/wine-stable/lib/wine/x86_64-unix",
    f"{WINE_DIR}/opt/wine-stable/lib/wine",
])

if len(sys.argv) < 2:
    print("usage: wine-appimage-run-pty.py <windows-program.exe> [args...]", file=sys.stderr)
    sys.exit(1)

env = os.environ.copy()
env.setdefault("WINEPREFIX", os.path.expanduser("~/.wine-appimage-stable"))
env["WINEDLLPATH"] = (
    f"{WINE_DIR}/opt/wine-stable/lib/wine/x86_64-unix:"
    f"{WINE_DIR}/opt/wine-stable/lib/wine"
)
env.setdefault("WINEDEBUG", "-all")

master, slave = pty.openpty()
proc = subprocess.Popen(
    [LOADER, "--library-path", LIBRARY_PATH,
     f"{WINE_DIR}/opt/wine-stable/bin/wine", *sys.argv[1:]],
    stdin=subprocess.PIPE, stdout=slave, stderr=slave, env=env,
    cwd=os.getcwd(), start_new_session=True,
)
os.close(slave)

def _terminate(_signum=None, _frame=None):
    try:
        os.killpg(proc.pid, signal.SIGTERM)
    except (ProcessLookupError, PermissionError):
        pass

signal.signal(signal.SIGTERM, _terminate)
signal.signal(signal.SIGINT, _terminate)

def _pump():
    answered = 0
    last_answer = 0.0
    tail = ""
    while True:
        r, _, _ = select.select([master], [], [], 1.0)
        if not r:
            if proc.poll() is not None:
                return
            continue
        try:
            data = os.read(master, 65536)
        except OSError:
            return
        if not data:
            return
        data = _ESC_SEQ.sub("", data.decode("utf-8", "replace"))
        if AUTO_ANSWER:
            tail = (tail + data)[-4096:]
            found = len(AUTO_PROMPT.findall(tail))
            if found > answered:
                answered = found
                now = time.monotonic()
                if now - last_answer >= 1.0:
                    last_answer = now
                    try:
                        proc.stdin.write((AUTO_ANSWER + "\n").encode())
                        proc.stdin.flush()
                        sys.stdout.write(f"\n[pty] auto-answer {AUTO_ANSWER!r}\n")
                        sys.stdout.flush()
                    except (BrokenPipeError, ValueError):
                        pass
        sys.stdout.write(data.replace("\r\n", "\n").replace("\r", "\n"))
        sys.stdout.flush()

try:
    _pump()
except KeyboardInterrupt:
    _terminate()

try:
    os.close(master)
except OSError:
    pass
try:
    proc.wait()
except Exception:
    pass
sys.exit(proc.returncode if proc.returncode is not None else 1)