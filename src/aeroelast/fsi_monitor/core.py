"""
Tolerant CSV reading for the FSI monitor.

The solver appends rows to ``rotor_performance.csv`` while the simulation
runs, so the file can be observed mid-write: trailing rows may be partial
and the file may be briefly missing. This module provides:

- :class:`SafeFileReader` — retrying, missing-file-tolerant raw file access.
- :class:`CSVReader` — cached, column-oriented access to a growing CSV file.

Both are dependency-free (stdlib ``csv`` only).
"""

from __future__ import annotations

import csv
import io
import threading
import time
from pathlib import Path
from typing import Any, List, Optional, Tuple

__all__ = ["CSVReader", "SafeFileReader"]


def _convert(value: str) -> Any:
    """Convert a CSV cell to ``float`` when possible, else keep the string."""
    v = value.strip()
    try:
        return float(v)
    except ValueError:
        return v


class SafeFileReader:
    """Tolerant reader for a file that may be actively written.

    Reads the raw file content, retrying transient I/O errors and returning
    ``None`` when the file is missing or unreadable. Instances hold no state
    about the file contents, so a single instance can be constructed once
    and reused for any number of files.

    Parameters
    ----------
    retries : int
        Number of retries for transient I/O errors before giving up.
    retry_delay : float
        Seconds to sleep between retries.
    encoding : str
        Encoding used to decode the file content.
    """

    def __init__(
        self,
        retries: int = 3,
        retry_delay: float = 0.05,
        encoding: str = "utf-8",
    ) -> None:
        self._retries = max(0, int(retries))
        self._retry_delay = max(0.0, retry_delay)
        self._encoding = encoding

    def read(self, path: Path) -> Optional[str]:
        """Return the full file content as text, or ``None`` on failure.

        A missing file returns ``None`` immediately. Transient OSErrors
        (e.g. a locked or half-flushed file) are retried up to ``retries``
        times before giving up.
        """
        target = Path(path)
        for attempt in range(self._retries + 1):
            try:
                return target.read_text(encoding=self._encoding)
            except FileNotFoundError:
                return None
            except OSError:
                if attempt < self._retries:
                    time.sleep(self._retry_delay)
                    continue
                return None
        return None


class CSVReader:
    """Cached, column-oriented reader for a growing CSV file.

    The file is re-read only when its size or modification time changes.
    Values are converted to ``float`` when possible; anything else is kept
    as ``str``. An unterminated final line (solver mid-write) is excluded
    until it is completed, and any row whose field count does not match the
    header is skipped without truncating the rest of the history.
    """

    def __init__(self, csv_path: Path, file_reader: SafeFileReader) -> None:
        self._csv_path = Path(csv_path)
        self._file_reader = file_reader

        self._lock = threading.Lock()
        self._columns: List[str] = []
        self._rows: List[dict[str, Any]] = []
        self._mtime_ns: Optional[int] = None
        self._size: Optional[int] = None

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    @property
    def columns(self) -> List[str]:
        """Header columns of the last successfully read file."""
        with self._lock:
            self._reload_if_changed()
            return list(self._columns)

    def get_latest_row(self) -> Optional[dict[str, Any]]:
        """Return the last complete data row as a ``dict``, or ``None``."""
        with self._lock:
            self._reload_if_changed()
            if not self._rows:
                return None
            return dict(self._rows[-1])

    def get_column(self, column: str, max_points: int = 0) -> List[Any]:
        """Return one column's values, tail-sliced when *max_points* > 0."""
        with self._lock:
            self._reload_if_changed()
            try:
                idx = self._columns.index(column)
            except ValueError:
                return []
            values = [row.get(self._columns[idx]) for row in self._rows]
            if max_points > 0:
                return values[-max_points:]
            return values

    def get_two_columns(
        self, x_col: str, y_col: str, max_points: int = 0
    ) -> Tuple[List[Any], List[Any]]:
        """Return paired (x, y) values, both tail-sliced consistently."""
        with self._lock:
            self._reload_if_changed()
            try:
                x_idx = self._columns.index(x_col)
                y_idx = self._columns.index(y_col)
            except ValueError:
                return [], []
            x_vals = [row.get(self._columns[x_idx]) for row in self._rows]
            y_vals = [row.get(self._columns[y_idx]) for row in self._rows]
            if max_points > 0:
                return x_vals[-max_points:], y_vals[-max_points:]
            return x_vals, y_vals

    # ------------------------------------------------------------------
    # Internal
    # ------------------------------------------------------------------

    def _reload_if_changed(self) -> None:
        """Re-parse the file when its size or mtime changed since last read."""
        try:
            stat = self._csv_path.stat()
        except OSError:
            return
        size = stat.st_size
        mtime_ns = stat.st_mtime_ns
        if size == self._size and mtime_ns == self._mtime_ns:
            return
        content = self._file_reader.read(self._csv_path)
        if content is None:
            return
        self._parse(content)
        self._size = size
        self._mtime_ns = mtime_ns

    def _parse(self, content: str) -> None:
        """Parse CSV text into header columns and data rows."""
        text = content
        if not content.endswith(("\n", "\r")):
            # Final line has no terminator -> the solver is likely mid-write.
            # Drop that line; it will be picked up once the write completes.
            newline_idx = content.rfind("\n")
            text = content[: newline_idx + 1] if newline_idx != -1 else ""
        reader = csv.reader(io.StringIO(text))
        columns: List[str] = []
        rows: List[dict[str, Any]] = []
        for raw in reader:
            if not raw:
                continue
            if not columns:
                columns = [col.strip() for col in raw]
                continue
            if len(raw) != len(columns):
                # Stale partial row left by a previous observation. Skip it
                # and keep going so later complete rows are not lost.
                continue
            rows.append({name: _convert(value) for name, value in zip(columns, raw)})
        self._columns = columns
        self._rows = rows
