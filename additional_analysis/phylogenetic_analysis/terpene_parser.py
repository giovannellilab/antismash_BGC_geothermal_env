from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import pandas as pd
from Bio import SeqIO
from Bio.SeqFeature import ExactPosition


# ============================================================
# PATHS — EDIT THESE
# ============================================================

BASE_DIR = Path(
    r"\..."
    r"\..."
)

MEMBERSHIP_CSV = BASE_DIR / "terpene_gcf_membership.csv"

MIBIG_GBK_ROOT = Path(
    r"\..."
    r"\..."
)

ENVIRONMENTAL_GBK_ROOT = Path(
    r"\..."
    r"\..."
)

# Main output used for MAFFT / trimming / IQ-TREE
COMBINED_FASTA = BASE_DIR / "terpene_tree_input_all_core_from_best_bgc.faa"

# Additional useful outputs
ENVIRONMENTAL_FASTA = (
    BASE_DIR / "terpene_environmental_best_bgc_all_core.faa"
)

MIBIG_FASTA = (
    BASE_DIR / "terpene_mibig_references.faa"
)

ENRICHED_OUTPUT = (
    BASE_DIR / "terpene_candidates_enriched.csv"
)

RANKED_ENVIRONMENTAL_OUTPUT = (
    BASE_DIR / "terpene_environmental_candidates_ranked.csv"
)

REPRESENTATIVES_OUTPUT = (
    BASE_DIR / "terpene_gcf_best_bgc_representatives.csv"
)

ENVIRONMENTAL_CORE_TABLE = (
    BASE_DIR / "terpene_environmental_best_bgc_all_core_proteins.csv"
)

MIBIG_CORE_TABLE = (
    BASE_DIR / "terpene_mibig_core_proteins.csv"
)

WARNINGS_OUTPUT = (
    BASE_DIR / "terpene_parser_warnings.txt"
)

# MIBiG reference BGCs that must be retained even when they are
# absent from terpene_gcf_membership.csv.
FORCED_MIBIG_IDS = {
    "BGC0002149",
    "BGC0000651",
    "BGC0001277",
}


# ============================================================
# CONTROLLED TERPENE CLASSIFICATION
# ============================================================

# Rules are evaluated from most specific to least specific.
# This avoids:
#   terpene_synth_c also matching terpene_synth
#   fung_ggpps2 also matching fung_ggpps
CLASSIFICATION_RULES: List[Tuple[str, List[str]]] = [
    ("fung_ggpps2", ["fung_ggpps2"]),
    ("terpene_synth_c", ["terpene_synth_c"]),
    ("trichodiene_synth", ["trichodiene_synth"]),
    ("terpene_cyclase", ["terpene_cyclase"]),
    ("phytoene_synth", ["phytoene_synt", "phytoene_synth"]),
    ("lycopene_cycl", ["lycopene_cycl"]),
    ("napt7", ["napt7"]),
    ("tri5", ["tri5"]),
    ("fung_ggpps", ["fung_ggpps"]),
    ("terpene_synth", ["terpene_synth"]),
]

# Used only to rank candidate BGCs, not to discard extra proteins.
# Lower number means higher priority.
ANNOTATION_PRIORITY = {
    "terpene_cyclase": 1,
    "terpene_synth": 1,
    "terpene_synth_c": 1,
    "trichodiene_synth": 1,
    "tri5": 1,

    "phytoene_synth": 2,
    "lycopene_cycl": 2,
    "napt7": 2,

    "fung_ggpps": 3,
    "fung_ggpps2": 3,

    "other_terpene_core": 50,
}


# ============================================================
# HELPERS
# ============================================================

def normalize_identifier(value: Any) -> str:
    """
    Convert text into a FASTA-safe identifier.
    """
    text = str(value).strip()

    for character in [
        " ",
        "|",
        ":",
        ";",
        ",",
        "(",
        ")",
        "[",
        "]",
        "{",
        "}",
        "/",
        "\\",
        "\t",
        "=",
    ]:
        text = text.replace(character, "_")

    while "__" in text:
        text = text.replace("__", "_")

    return text.strip("_")


def clean_translation(translation: Any) -> str:
    """
    Remove whitespace and a terminal stop symbol.
    """
    sequence = (
        str(translation)
        .replace(" ", "")
        .replace("\n", "")
        .replace("\r", "")
        .replace("\t", "")
    )

    return sequence.rstrip("*")


def parse_boolean(value: Optional[str]) -> Optional[bool]:
    if value is None:
        return None

    normalized = str(value).strip().lower()

    if normalized == "true":
        return True

    if normalized == "false":
        return False

    return None


def qualifier_text(feature: Any) -> str:
    """
    Combine antiSMASH annotation fields into one searchable string.
    """
    gene_functions = " ".join(
        feature.qualifiers.get("gene_functions", [])
    )

    sec_met_domain = " ".join(
        feature.qualifiers.get("sec_met_domain", [])
    )

    gene = " ".join(
        feature.qualifiers.get("gene", [])
    )

    product = " ".join(
        feature.qualifiers.get("product", [])
    )

    return (
        f"{gene_functions} {sec_met_domain} {gene} {product}"
        .lower()
    )


def classify_feature(feature: Any) -> Optional[str]:
    """
    Return one controlled terpene classification.

    Returns None when none of the requested terpene keywords are present.
    """
    text = qualifier_text(feature)

    for classification, patterns in CLASSIFICATION_RULES:
        if any(pattern in text for pattern in patterns):
            return classification

    return None


def annotation_priority(classification: Optional[str]) -> int:
    if classification is None:
        return 99

    return ANNOTATION_PRIORITY.get(
        classification,
        ANNOTATION_PRIORITY["other_terpene_core"],
    )


def feature_touches_record_edge(
    feature: Any,
    record_length: int,
) -> bool:
    start = int(feature.location.start)
    end = int(feature.location.end)

    return start == 0 or end == record_length


def feature_has_fuzzy_location(feature: Any) -> bool:
    return not (
        isinstance(feature.location.start, ExactPosition)
        and isinstance(feature.location.end, ExactPosition)
    )


def find_region_feature(record: Any) -> Optional[Any]:
    for feature in record.features:
        if feature.type == "region":
            return feature

    return None


def safe_bool_column(
    dataframe: pd.DataFrame,
    column: str,
    default: bool,
) -> pd.Series:
    if column not in dataframe.columns:
        return pd.Series(
            [default] * len(dataframe),
            index=dataframe.index,
            dtype=bool,
        )

    return dataframe[column].fillna(default).astype(bool)


def safe_numeric_column(
    dataframe: pd.DataFrame,
    column: str,
    default: float = 0,
) -> pd.Series:
    if column not in dataframe.columns:
        return pd.Series(
            [default] * len(dataframe),
            index=dataframe.index,
            dtype=float,
        )

    return pd.to_numeric(
        dataframe[column],
        errors="coerce",
    ).fillna(default)


def write_wrapped_sequence(
    handle: Any,
    sequence: str,
    width: int = 80,
) -> None:
    for start in range(0, len(sequence), width):
        handle.write(sequence[start:start + width] + "\n")


# ============================================================
# GBK FILE LOCATION
# ============================================================

def build_filename_index(
    root: Path,
) -> Dict[str, List[Path]]:
    index: Dict[str, List[Path]] = {}

    for path in root.rglob("*.gbk"):
        index.setdefault(path.name, []).append(path)

    return index


def resolve_local_gbk(
    filename: str,
    filename_index: Dict[str, List[Path]],
    expected_path: Optional[str] = None,
) -> Tuple[Optional[Path], Optional[str]]:
    matches = filename_index.get(filename, [])

    if not matches:
        return None, f"GBK not found: {filename}"

    if len(matches) == 1:
        return matches[0], None

    if expected_path:
        expected_normalized = (
            str(expected_path)
            .replace("\\", "/")
            .lower()
        )

        scored_matches: List[Tuple[int, Path]] = []

        for path in matches:
            local_normalized = (
                str(path)
                .replace("\\", "/")
                .lower()
            )

            score = 0

            for component in Path(expected_normalized).parts:
                component_text = str(component).lower()

                if (
                    component_text
                    and component_text in local_normalized
                ):
                    score += 1

            scored_matches.append((score, path))

        scored_matches.sort(
            key=lambda item: (-item[0], str(item[1]))
        )

        best_score, best_path = scored_matches[0]

        if best_score > 0:
            return (
                best_path,
                (
                    f"Multiple GBKs named {filename}; "
                    f"selected {best_path} using path similarity."
                ),
            )

    all_paths = " | ".join(str(path) for path in matches)

    return (
        matches[0],
        (
            f"Multiple GBKs named {filename}; "
            f"using first match: {matches[0]}. "
            f"All matches: {all_paths}"
        ),
    )


def resolve_forced_mibig_gbk(
    mibig_id: str,
    filename_index: Dict[str, List[Path]],
) -> Tuple[Optional[Path], Optional[str]]:
    """
    Locate a forced MIBiG GBK by accession.

    Exact filenames such as BGC0002149.gbk are preferred. Files whose
    names begin with the accession are accepted as a fallback.
    """
    exact_filename = f"{mibig_id}.gbk"
    exact_matches = filename_index.get(exact_filename, [])

    if exact_matches:
        if len(exact_matches) == 1:
            return exact_matches[0], None

        return (
            sorted(exact_matches, key=str)[0],
            (
                f"Multiple exact GBKs found for forced reference "
                f"{mibig_id}; using "
                f"{sorted(exact_matches, key=str)[0]}."
            ),
        )

    prefix_matches: List[Path] = []

    for filename, paths in filename_index.items():
        if filename.upper().startswith(mibig_id.upper()):
            prefix_matches.extend(paths)

    prefix_matches = sorted(set(prefix_matches), key=str)

    if not prefix_matches:
        return (
            None,
            f"Forced MIBiG GBK not found for accession {mibig_id}.",
        )

    if len(prefix_matches) == 1:
        return prefix_matches[0], None

    return (
        prefix_matches[0],
        (
            f"Multiple prefix-matching GBKs found for forced reference "
            f"{mibig_id}; using {prefix_matches[0]}. "
            f"All matches: {' | '.join(str(path) for path in prefix_matches)}"
        ),
    )


# ============================================================
# GBK PARSING
# ============================================================

def parse_gbk(
    gbk_path: Path,
) -> Dict[str, Any]:
    """
    Parse all records from one antiSMASH GBK.

    All matching terpene-associated CDSs are retained.
    A single best CDS is also recorded only for BGC-ranking purposes.
    """
    records = list(SeqIO.parse(str(gbk_path), "genbank"))

    if not records:
        raise ValueError("No GenBank records found")

    all_core_candidates: List[Dict[str, Any]] = []
    total_n_cds = 0
    total_length = 0

    record_ids: List[str] = []
    descriptions: List[str] = []
    region_products: List[str] = []
    region_numbers: List[str] = []
    contig_edge_values: List[bool] = []

    for record_index, record in enumerate(records, start=1):
        record_length = len(record.seq)
        total_length += record_length
        record_ids.append(str(record.id))
        descriptions.append(str(record.description))

        region_feature = find_region_feature(record)

        if region_feature is not None:
            edge_values = region_feature.qualifiers.get(
                "contig_edge",
                [None],
            )

            edge_value = edge_values[0] if edge_values else None
            parsed_edge = parse_boolean(edge_value)

            if parsed_edge is not None:
                contig_edge_values.append(parsed_edge)

            region_products.extend(
                region_feature.qualifiers.get("product", [])
            )

            region_numbers.extend(
                region_feature.qualifiers.get("region_number", [])
            )

        for feature_index, feature in enumerate(
            record.features,
            start=1,
        ):
            if feature.type != "CDS":
                continue

            total_n_cds += 1

            classification = classify_feature(feature)

            if classification is None:
                continue

            translations = feature.qualifiers.get(
                "translation",
                [],
            )

            if not translations:
                continue

            translation = clean_translation(translations[0])

            if not translation:
                continue

            locus_tag = feature.qualifiers.get(
                "locus_tag",
                ["unknown"],
            )[0]

            protein_id = feature.qualifiers.get(
                "protein_id",
                [""],
            )[0]

            gene_name = feature.qualifiers.get(
                "gene",
                [""],
            )[0]

            gene_functions = " | ".join(
                feature.qualifiers.get(
                    "gene_functions",
                    [],
                )
            )

            sec_met_domain = " | ".join(
                feature.qualifiers.get(
                    "sec_met_domain",
                    [],
                )
            )

            gene_kind = " | ".join(
                feature.qualifiers.get(
                    "gene_kind",
                    [],
                )
            )

            product = " | ".join(
                feature.qualifiers.get(
                    "product",
                    [],
                )
            )

            touches_edge = feature_touches_record_edge(
                feature,
                record_length,
            )

            fuzzy_location = feature_has_fuzzy_location(
                feature
            )

            likely_non_truncated = (
                not touches_edge
                and not fuzzy_location
                and len(translation) >= 100
            )

            candidate = {
                "candidate_index": len(all_core_candidates) + 1,
                "record_index": record_index,
                "record_id": str(record.id),
                "feature_index": feature_index,
                "locus_tag": str(locus_tag),
                "protein_id": str(protein_id),
                "gene_name": str(gene_name),
                "classification": classification,
                "classification_priority": annotation_priority(
                    classification
                ),
                "translation": translation,
                "core_length_aa": len(translation),
                "gene_functions": gene_functions,
                "sec_met_domain": sec_met_domain,
                "gene_kind": gene_kind,
                "product": product,
                "core_start": int(feature.location.start),
                "core_end": int(feature.location.end),
                "core_strand": feature.location.strand,
                "core_touches_record_edge": touches_edge,
                "core_fuzzy_location": fuzzy_location,
                "core_likely_non_truncated": likely_non_truncated,
            }

            all_core_candidates.append(candidate)

    ranked_core_candidates = sorted(
        all_core_candidates,
        key=lambda candidate: (
            not candidate["core_likely_non_truncated"],
            candidate["classification_priority"],
            -candidate["core_length_aa"],
            candidate["locus_tag"],
            candidate["candidate_index"],
        ),
    )

    best_core = (
        ranked_core_candidates[0]
        if ranked_core_candidates
        else None
    )

    # Conservative BGC-level edge status:
    # True if any parsed record says it is at a contig edge.
    contig_edge: Optional[bool]

    if contig_edge_values:
        contig_edge = any(contig_edge_values)
    else:
        contig_edge = None

    result: Dict[str, Any] = {
        "parsed_gbk_path": str(gbk_path),
        "record_id": ";".join(record_ids),
        "record_description": " | ".join(descriptions),
        "bgc_length_bp": total_length,
        "n_cds": total_n_cds,
        "contig_edge": contig_edge,
        "region_product_from_gbk": ";".join(region_products),
        "region_number": ";".join(region_numbers),
        "n_core_candidates": len(all_core_candidates),
        "has_valid_core": bool(all_core_candidates),
        "multiple_core_candidates": len(all_core_candidates) > 1,
        "all_core_candidates": all_core_candidates,

        # Best CDS is used only to rank BGCs.
        "selected_core_locus": None,
        "selected_core_protein_id": None,
        "selected_core_gene_name": None,
        "selected_core_classification": None,
        "selected_core_gene_functions": None,
        "selected_core_sec_met_domain": None,
        "selected_core_gene_kind": None,
        "selected_core_start": None,
        "selected_core_end": None,
        "selected_core_strand": None,
        "core_length_aa": None,
        "core_touches_record_edge": None,
        "core_fuzzy_location": None,
        "core_likely_non_truncated": False,
        "selected_translation": None,
    }

    if best_core is not None:
        result.update(
            {
                "selected_core_locus": best_core["locus_tag"],
                "selected_core_protein_id": best_core["protein_id"],
                "selected_core_gene_name": best_core["gene_name"],
                "selected_core_classification": best_core[
                    "classification"
                ],
                "selected_core_gene_functions": best_core[
                    "gene_functions"
                ],
                "selected_core_sec_met_domain": best_core[
                    "sec_met_domain"
                ],
                "selected_core_gene_kind": best_core[
                    "gene_kind"
                ],
                "selected_core_start": best_core["core_start"],
                "selected_core_end": best_core["core_end"],
                "selected_core_strand": best_core["core_strand"],
                "core_length_aa": best_core["core_length_aa"],
                "core_touches_record_edge": best_core[
                    "core_touches_record_edge"
                ],
                "core_fuzzy_location": best_core[
                    "core_fuzzy_location"
                ],
                "core_likely_non_truncated": best_core[
                    "core_likely_non_truncated"
                ],
                "selected_translation": best_core["translation"],
            }
        )

    return result


# ============================================================
# MAIN WORKFLOW
# ============================================================

def main() -> None:
    if not MEMBERSHIP_CSV.exists():
        raise FileNotFoundError(
            f"Membership CSV not found:\n{MEMBERSHIP_CSV}"
        )

    if not MIBIG_GBK_ROOT.exists():
        raise FileNotFoundError(
            f"MIBiG GBK directory not found:\n{MIBIG_GBK_ROOT}"
        )

    if not ENVIRONMENTAL_GBK_ROOT.exists():
        raise FileNotFoundError(
            "Environmental GBK directory not found:\n"
            f"{ENVIRONMENTAL_GBK_ROOT}"
        )

    BASE_DIR.mkdir(parents=True, exist_ok=True)

    print("Reading BGC-to-GCF membership table...")

    membership = pd.read_csv(MEMBERSHIP_CSV)

    required_columns = {
        "gcf_id",
        "bgc_record_id",
        "product",
        "source_type",
        "biosample_or_mibig_id",
        "gbk_filename",
        "gbk_path",
    }

    missing_columns = required_columns.difference(
        membership.columns
    )

    if missing_columns:
        raise ValueError(
            "Membership table is missing required columns: "
            + ", ".join(sorted(missing_columns))
        )

    for column in [
        "source_type",
        "product",
        "gbk_filename",
        "biosample_or_mibig_id",
    ]:
        membership[column] = (
            membership[column]
            .astype(str)
            .str.strip()
        )

    print("Indexing environmental GBKs...")
    environmental_filename_index = build_filename_index(
        ENVIRONMENTAL_GBK_ROOT
    )

    print("Indexing MIBiG GBKs...")
    mibig_filename_index = build_filename_index(
        MIBIG_GBK_ROOT
    )

    warnings: List[str] = []
    parsed_rows: List[Dict[str, Any]] = []

    # Stores all CDS candidates for each parsed membership row.
    candidate_map: Dict[
        Tuple[str, str, str],
        List[Dict[str, Any]],
    ] = {}

    total_rows = len(membership)

    print("Parsing membership GBKs...")

    for row_number, row in enumerate(
        membership.itertuples(index=False),
        start=1,
    ):
        base_data = row._asdict()

        if row_number % 100 == 0 or row_number == total_rows:
            print(f"Processed {row_number}/{total_rows}")

        source_type = str(row.source_type).strip().lower()
        expected_path = str(getattr(row, "gbk_path", ""))

        if source_type == "mibig":
            filename_index = mibig_filename_index
        elif source_type == "environmental":
            filename_index = environmental_filename_index
        else:
            parsed_rows.append(
                {
                    **base_data,
                    "file_found": False,
                    "parse_status": "unknown_source_type",
                }
            )

            warnings.append(
                f"Unknown source_type '{row.source_type}' "
                f"for {row.gbk_filename}"
            )
            continue

        local_path, warning = resolve_local_gbk(
            filename=str(row.gbk_filename),
            filename_index=filename_index,
            expected_path=expected_path,
        )

        if warning:
            warnings.append(warning)

        if local_path is None:
            parsed_rows.append(
                {
                    **base_data,
                    "file_found": False,
                    "parse_status": "file_not_found",
                }
            )
            continue

        try:
            parsed_data = parse_gbk(local_path)

            all_core_candidates = parsed_data.pop(
                "all_core_candidates",
                [],
            )

            row_key = (
                str(row.gcf_id),
                str(row.bgc_record_id),
                str(local_path),
            )

            candidate_map[row_key] = all_core_candidates

            parsed_rows.append(
                {
                    **base_data,
                    "file_found": True,
                    "parse_status": "parsed",
                    **parsed_data,
                }
            )

        except Exception as error:
            warnings.append(
                f"Could not parse {local_path}: {error}"
            )

            parsed_rows.append(
                {
                    **base_data,
                    "file_found": True,
                    "parse_status": f"parse_error: {error}",
                }
            )

    enriched = pd.DataFrame(parsed_rows)

    if enriched.empty:
        raise ValueError("No membership records were processed.")

    enriched.to_csv(ENRICHED_OUTPUT, index=False)

    # ========================================================
    # SELECT ONE BEST ENVIRONMENTAL BGC PER GCF
    # ========================================================

    environmental = enriched[
        enriched["source_type"]
        .astype(str)
        .str.strip()
        .str.lower()
        .eq("environmental")
    ].copy()

    if environmental.empty:
        raise ValueError(
            "No environmental rows were found."
        )

    environmental["has_valid_core"] = safe_bool_column(
        environmental,
        "has_valid_core",
        False,
    )

    environmental["contig_edge"] = safe_bool_column(
        environmental,
        "contig_edge",
        True,
    )

    environmental[
        "core_likely_non_truncated"
    ] = safe_bool_column(
        environmental,
        "core_likely_non_truncated",
        False,
    )

    environmental["core_length_aa"] = safe_numeric_column(
        environmental,
        "core_length_aa",
        0,
    )

    environmental["n_core_candidates"] = safe_numeric_column(
        environmental,
        "n_core_candidates",
        0,
    )

    environmental["n_cds"] = safe_numeric_column(
        environmental,
        "n_cds",
        0,
    )

    environmental["bgc_length_bp"] = safe_numeric_column(
        environmental,
        "bgc_length_bp",
        0,
    )

    environmental["bgc_record_id"] = pd.to_numeric(
        environmental["bgc_record_id"],
        errors="raise",
    )

    environmental["product_priority"] = (
        environmental["product"]
        .astype(str)
        .str.strip()
        .str.lower()
        .map(
            {
                "terpene": 1,
                "terpene-precursor": 2,
            }
        )
        .fillna(3)
    )

    valid_environmental = environmental[
        environmental["parse_status"].astype(str).eq("parsed")
        &
        environmental["file_found"].fillna(False).astype(bool)
        &
        environmental["has_valid_core"]
    ].copy()

    # Ranking order:
    # 1. complete BGCs before contig-edge BGCs
    # 2. non-truncated best core protein
    # 3. terpene product class
    # 4. more terpene-core proteins in the BGC
    # 5. longer best core protein
    # 6. more CDSs
    # 7. longer BGC
    valid_environmental = valid_environmental.sort_values(
        by=[
            "gcf_id",
            "contig_edge",
            "core_likely_non_truncated",
            "product_priority",
            "n_core_candidates",
            "core_length_aa",
            "n_cds",
            "bgc_length_bp",
            "bgc_record_id",
        ],
        ascending=[
            True,
            True,
            False,
            True,
            False,
            False,
            False,
            False,
            True,
        ],
    )

    valid_environmental["representative_rank"] = (
        valid_environmental
        .groupby("gcf_id")
        .cumcount()
        + 1
    )

    valid_environmental.to_csv(
        RANKED_ENVIRONMENTAL_OUTPUT,
        index=False,
    )

    representatives = valid_environmental[
        valid_environmental["representative_rank"] == 1
    ].copy()

    representatives.to_csv(
        REPRESENTATIVES_OUTPUT,
        index=False,
    )

    # ========================================================
    # WRITE ALL CORE PROTEINS FROM EACH SELECTED BGC
    # ========================================================

    environmental_fasta_records: List[
        Tuple[str, str]
    ] = []

    environmental_core_rows: List[
        Dict[str, Any]
    ] = []

    seen_environmental_headers: Dict[str, int] = {}

    for row in representatives.itertuples(index=False):
        row_key = (
            str(row.gcf_id),
            str(row.bgc_record_id),
            str(row.parsed_gbk_path),
        )

        candidates = candidate_map.get(row_key, [])

        if not candidates:
            warnings.append(
                f"GCF {row.gcf_id}: selected BGC has no "
                "recoverable candidate list."
            )
            continue

        biosample = normalize_identifier(
            row.biosample_or_mibig_id
        )

        for candidate in candidates:
            sequence = clean_translation(
                candidate["translation"]
            )

            if not sequence:
                continue

            locus = normalize_identifier(
                candidate.get("locus_tag", "unknown")
            )

            classification = normalize_identifier(
                candidate.get(
                    "classification",
                    "other_terpene_core",
                )
            )

            base_header = (
                f"ENV_GCF_{row.gcf_id}"
                f"__{biosample}"
                f"__REC_{row.bgc_record_id}"
                f"__{locus}"
                f"__{classification}"
            )

            seen_environmental_headers[base_header] = (
                seen_environmental_headers.get(
                    base_header,
                    0,
                )
                + 1
            )

            occurrence = seen_environmental_headers[
                base_header
            ]

            tip_id = (
                base_header
                if occurrence == 1
                else f"{base_header}__DUP_{occurrence}"
            )

            environmental_fasta_records.append(
                (tip_id, sequence)
            )

            environmental_core_rows.append(
                {
                    "tip_id": tip_id,
                    "gcf_id": row.gcf_id,
                    "bgc_record_id": row.bgc_record_id,
                    "biosample_or_mibig_id": (
                        row.biosample_or_mibig_id
                    ),
                    "gbk_filename": row.gbk_filename,
                    "parsed_gbk_path": row.parsed_gbk_path,
                    **{
                        key: value
                        for key, value in candidate.items()
                        if key != "translation"
                    },
                    "translation": sequence,
                }
            )

    pd.DataFrame(
        environmental_core_rows
    ).to_csv(
        ENVIRONMENTAL_CORE_TABLE,
        index=False,
    )

    with ENVIRONMENTAL_FASTA.open(
        "w",
        encoding="utf-8",
    ) as handle:
        for tip_id, sequence in environmental_fasta_records:
            handle.write(f">{tip_id}\n")
            write_wrapped_sequence(handle, sequence)

    # ========================================================
    # WRITE ALL MIBIG REFERENCE CORE PROTEINS
    # ========================================================

    mibig_rows = enriched[
        enriched["source_type"]
        .astype(str)
        .str.strip()
        .str.lower()
        .eq("mibig")
        &
        enriched["parse_status"]
        .astype(str)
        .eq("parsed")
    ].copy()

    mibig_fasta_records: List[
        Tuple[str, str]
    ] = []

    mibig_core_rows: List[
        Dict[str, Any]
    ] = []

    seen_mibig_headers: Dict[str, int] = {}
    included_mibig_ids: set[str] = set()

    for row in mibig_rows.itertuples(index=False):
        row_key = (
            str(row.gcf_id),
            str(row.bgc_record_id),
            str(row.parsed_gbk_path),
        )

        candidates = candidate_map.get(row_key, [])

        mibig_id = normalize_identifier(
            row.biosample_or_mibig_id
        )

        for candidate in candidates:
            sequence = clean_translation(
                candidate["translation"]
            )

            if not sequence:
                continue

            locus = normalize_identifier(
                candidate.get("locus_tag", "unknown")
            )

            classification = normalize_identifier(
                candidate.get(
                    "classification",
                    "other_terpene_core",
                )
            )

            base_header = (
                f"MIBIG_{mibig_id}"
                f"__GCF_{row.gcf_id}"
                f"__{locus}"
                f"__{classification}"
            )

            seen_mibig_headers[base_header] = (
                seen_mibig_headers.get(base_header, 0)
                + 1
            )

            occurrence = seen_mibig_headers[base_header]

            tip_id = (
                base_header
                if occurrence == 1
                else f"{base_header}__DUP_{occurrence}"
            )

            mibig_fasta_records.append(
                (tip_id, sequence)
            )

            included_mibig_ids.add(
                str(row.biosample_or_mibig_id).strip().upper()
            )

            mibig_core_rows.append(
                {
                    "tip_id": tip_id,
                    "gcf_id": row.gcf_id,
                    "bgc_record_id": row.bgc_record_id,
                    "biosample_or_mibig_id": (
                        row.biosample_or_mibig_id
                    ),
                    "gbk_filename": row.gbk_filename,
                    "parsed_gbk_path": row.parsed_gbk_path,
                    **{
                        key: value
                        for key, value in candidate.items()
                        if key != "translation"
                    },
                    "translation": sequence,
                }
            )

    # ========================================================
    # FORCE-RETAIN SPECIFIED MIBIG REFERENCE BGCs
    # ========================================================

    missing_forced_references: List[str] = []

    for forced_mibig_id in sorted(FORCED_MIBIG_IDS):
        forced_mibig_id_upper = forced_mibig_id.upper()

        # Do not duplicate a reference already exported normally.
        if forced_mibig_id_upper in included_mibig_ids:
            continue

        forced_path, warning = resolve_forced_mibig_gbk(
            mibig_id=forced_mibig_id,
            filename_index=mibig_filename_index,
        )

        if warning:
            warnings.append(warning)

        if forced_path is None:
            missing_forced_references.append(
                f"{forced_mibig_id}: GBK file not found"
            )
            continue

        try:
            forced_parsed_data = parse_gbk(forced_path)
        except Exception as error:
            warnings.append(
                f"Could not parse forced reference "
                f"{forced_mibig_id} ({forced_path}): {error}"
            )
            missing_forced_references.append(
                f"{forced_mibig_id}: parse error"
            )
            continue

        forced_candidates = forced_parsed_data.get(
            "all_core_candidates",
            [],
        )

        if not forced_candidates:
            warnings.append(
                f"Forced reference {forced_mibig_id} was found and "
                "parsed, but none of its CDS features matched the "
                "controlled terpene classification rules."
            )
            missing_forced_references.append(
                f"{forced_mibig_id}: no matching terpene-core CDS"
            )
            continue

        for candidate in forced_candidates:
            sequence = clean_translation(
                candidate["translation"]
            )

            if not sequence:
                continue

            locus = normalize_identifier(
                candidate.get("locus_tag", "unknown")
            )

            classification = normalize_identifier(
                candidate.get(
                    "classification",
                    "other_terpene_core",
                )
            )

            base_header = (
                f"MIBIG_{forced_mibig_id}"
                f"__GCF_unassigned"
                f"__{locus}"
                f"__{classification}"
            )

            seen_mibig_headers[base_header] = (
                seen_mibig_headers.get(base_header, 0)
                + 1
            )

            occurrence = seen_mibig_headers[base_header]

            tip_id = (
                base_header
                if occurrence == 1
                else f"{base_header}__DUP_{occurrence}"
            )

            mibig_fasta_records.append(
                (tip_id, sequence)
            )

            included_mibig_ids.add(forced_mibig_id_upper)

            mibig_core_rows.append(
                {
                    "tip_id": tip_id,
                    "gcf_id": "unassigned",
                    "bgc_record_id": "",
                    "biosample_or_mibig_id": forced_mibig_id,
                    "gbk_filename": forced_path.name,
                    "parsed_gbk_path": str(forced_path),
                    "forced_reference": True,
                    **{
                        key: value
                        for key, value in candidate.items()
                        if key != "translation"
                    },
                    "translation": sequence,
                }
            )

    # Stop rather than silently creating another tree without one of
    # the three reviewer-relevant anchor references.
    if missing_forced_references:
        WARNINGS_OUTPUT.write_text(
            "\n".join(warnings),
            encoding="utf-8",
        )

        raise RuntimeError(
            "The following forced MIBiG references could not be "
            "included:\n- "
            + "\n- ".join(missing_forced_references)
            + "\n\nCheck the GBK filenames and, when a GBK was "
            "parsed but had no matching CDS, inspect its annotations "
            "and expand CLASSIFICATION_RULES before rebuilding the tree."
        )

    pd.DataFrame(
        mibig_core_rows
    ).to_csv(
        MIBIG_CORE_TABLE,
        index=False,
    )

    with MIBIG_FASTA.open(
        "w",
        encoding="utf-8",
    ) as handle:
        for tip_id, sequence in mibig_fasta_records:
            handle.write(f">{tip_id}\n")
            write_wrapped_sequence(handle, sequence)

    # ========================================================
    # COMBINED FASTA FOR ALIGNMENT AND IQ-TREE
    # ========================================================

    with COMBINED_FASTA.open(
        "w",
        encoding="utf-8",
    ) as handle:
        for tip_id, sequence in environmental_fasta_records:
            handle.write(f">{tip_id}\n")
            write_wrapped_sequence(handle, sequence)

        for tip_id, sequence in mibig_fasta_records:
            handle.write(f">{tip_id}\n")
            write_wrapped_sequence(handle, sequence)

    # ========================================================
    # WARNINGS AND SUMMARY
    # ========================================================

    all_environmental_gcfs = set(
        environmental["gcf_id"]
        .dropna()
        .unique()
    )

    selected_environmental_gcfs = set(
        representatives["gcf_id"]
        .dropna()
        .unique()
    )

    missing_environmental_gcfs = sorted(
        all_environmental_gcfs
        - selected_environmental_gcfs
    )

    for gcf_id in missing_environmental_gcfs:
        warnings.append(
            f"GCF {gcf_id}: no successfully parsed "
            "environmental BGC with a valid terpene core."
        )

    WARNINGS_OUTPUT.write_text(
        "\n".join(warnings),
        encoding="utf-8",
    )

    print()
    print("=" * 72)
    print("FINISHED")
    print("=" * 72)
    print(f"Membership rows: {len(membership)}")
    print(
        "Environmental representative BGCs selected: "
        f"{len(representatives)}"
    )
    print(
        "Environmental core proteins retained from selected BGCs: "
        f"{len(environmental_fasta_records)}"
    )
    print(
        "MIBiG reference proteins retained: "
        f"{len(mibig_fasta_records)}"
    )
    print(
        "Forced MIBiG references confirmed: "
        + ", ".join(sorted(FORCED_MIBIG_IDS))
    )
    print(
        "Combined sequences for tree: "
        f"{len(environmental_fasta_records) + len(mibig_fasta_records)}"
    )
    print(
        "Environmental GCFs without a valid representative: "
        f"{len(missing_environmental_gcfs)}"
    )
    print(f"Warnings: {len(warnings)}")
    print()
    print(f"Combined FASTA:\n{COMBINED_FASTA}")
    print(f"Warnings:\n{WARNINGS_OUTPUT}")


if __name__ == "__main__":
    main()