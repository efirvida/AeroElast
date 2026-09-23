"""Helpers for consuming generated Aero-FSI report manifests."""

from __future__ import annotations

import json
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from types import MappingProxyType
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import polars as pl

DEFAULT_AERO_REPORT_SCHEMA_MANIFEST_FILENAME = "aero_report_schema.json"


def _require_mapping(value: object, *, context: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise ValueError(f"{context} must be a mapping")
    return value


def _require_str(mapping: Mapping[str, object], key: str, *, context: str) -> str:
    value = mapping.get(key)
    if not isinstance(value, str) or not value:
        raise ValueError(f"{context}.{key} must be a non-empty string")
    return value


def _require_int(mapping: Mapping[str, object], key: str, *, context: str) -> int:
    value = mapping.get(key)
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError(f"{context}.{key} must be an integer")
    return value


def _require_fieldnames(
    mapping: Mapping[str, object], key: str, *, context: str
) -> tuple[str, ...]:
    value = mapping.get(key)
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)):
        raise ValueError(f"{context}.{key} must be a sequence of strings")

    fieldnames = tuple(str(field) for field in value)
    if not fieldnames or any(not field for field in fieldnames):
        raise ValueError(f"{context}.{key} must contain non-empty strings")
    return fieldnames


@dataclass(frozen=True, slots=True)
class AeroReportDefinition:
    report_type: str
    path: Path
    schema_version: str
    granularity: str
    summary: str
    fieldnames: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class AeroReportSchemaManifest:
    output_folder: Path
    manifest_path: Path
    manifest_version: str
    backend: str
    n_blades: int
    sectional_bins: int
    reports: Mapping[str, AeroReportDefinition]

    def get_report(self, report_type: str) -> AeroReportDefinition:
        try:
            return self.reports[report_type]
        except KeyError as exc:
            raise KeyError(f"Unknown Aero-FSI report type: {report_type}") from exc

    @property
    def aero_report(self) -> AeroReportDefinition:
        return self.get_report("aero_report")

    @property
    def aero_spanwise_report(self) -> AeroReportDefinition:
        return self.get_report("aero_spanwise_report")

    @property
    def aero_sectional_report(self) -> AeroReportDefinition:
        return self.get_report("aero_sectional_report")


def load_aero_report_schema(
    output_folder: str | Path,
    *,
    manifest_filename: str = DEFAULT_AERO_REPORT_SCHEMA_MANIFEST_FILENAME,
) -> AeroReportSchemaManifest:
    """Load the generated Aero-FSI report schema manifest from an output folder."""
    output_path = Path(output_folder)
    manifest_path = output_path / manifest_filename
    if not manifest_path.exists():
        raise FileNotFoundError(f"Aero-FSI report schema manifest not found: {manifest_path}")

    with manifest_path.open(encoding="utf-8") as handle:
        payload = json.load(handle)

    manifest = _require_mapping(payload, context="manifest")
    reports_payload = _require_mapping(manifest.get("reports"), context="manifest.reports")

    report_map: dict[str, AeroReportDefinition] = {}
    for report_type, report_payload in reports_payload.items():
        report_context = f"manifest.reports.{report_type}"
        report = _require_mapping(report_payload, context=report_context)
        report_file = Path(_require_str(report, "file", context=report_context))
        report_map[str(report_type)] = AeroReportDefinition(
            report_type=str(report_type),
            path=report_file if report_file.is_absolute() else output_path / report_file,
            schema_version=_require_str(report, "schema_version", context=report_context),
            granularity=_require_str(report, "granularity", context=report_context),
            summary=_require_str(report, "summary", context=report_context),
            fieldnames=_require_fieldnames(report, "fieldnames", context=report_context),
        )

    return AeroReportSchemaManifest(
        output_folder=output_path,
        manifest_path=manifest_path,
        manifest_version=_require_str(manifest, "manifest_version", context="manifest"),
        backend=_require_str(manifest, "backend", context="manifest"),
        n_blades=_require_int(manifest, "n_blades", context="manifest"),
        sectional_bins=_require_int(manifest, "sectional_bins", context="manifest"),
        reports=MappingProxyType(report_map),
    )


def load_aero_report_table(
    output_folder: str | Path,
    report_type: str,
    *,
    validate_columns: bool = True,
) -> "pl.DataFrame":
    """Load one generated Aero-FSI CSV and validate its columns against the manifest."""
    import polars as pl

    schema = load_aero_report_schema(output_folder)
    report = schema.get_report(report_type)
    if not report.path.exists():
        raise FileNotFoundError(f"Aero-FSI report file not found: {report.path}")

    table = pl.read_csv(report.path)
    if validate_columns and tuple(table.columns) != report.fieldnames:
        raise ValueError(
            "Aero-FSI report columns do not match the manifest field order for "
            f"{report_type}: {report.path}"
        )
    return table