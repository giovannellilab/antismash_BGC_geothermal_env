from pathlib import Path
import pandas as pd
from pathlib import Path

THRESHOLD = 0.70  # keep only contigs with genus score > this

def build_genus_dict_from_tsv(tsv_path: str | Path, threshold: float = THRESHOLD) -> dict:
    """
    Reads your tab-delimited classification file and returns:
      {contig_id: {"genus": <name>, "genus_percentage": <float>}}
    Only includes entries with genus_percentage > threshold.
    Robust to ragged rows (e.g. 'no taxid assigned' lines).
    """
    df = pd.read_csv(
        tsv_path,
        sep="\t",
        comment="#",
        header=None,
        names=["contig", "classification", "reason", "lineage", "scores"],
        engine="python"  # important for ragged rows
    )

    result = {}

    for _, row in df.iterrows():
        contig = row.get("contig")
        lineage = row.get("lineage")
        scores = row.get("scores")

        if pd.isna(contig) or pd.isna(lineage) or pd.isna(scores):
            continue

        lineage_parts = str(lineage).split(";")

        # parse score list safely
        try:
            score_parts = [float(x) for x in str(scores).split(";") if x != ""]
        except ValueError:
            continue

        genus_idx = None
        genus_name = None

        for i, part in enumerate(lineage_parts):
            if part.startswith("g__"):
                genus_idx = i
                genus_name = part.replace("g__", "")
                break

        if genus_idx is None or genus_idx >= len(score_parts):
            continue

        genus_score = score_parts[genus_idx]

        if genus_score > threshold:
            result[str(contig)] = {
                "genus": genus_name,
                "genus_percentage": genus_score
            }

    return result


def add_contig_lengths_from_fasta(fasta_path: str | Path, contig_dict: dict) -> dict:
    """
    For each contig in the FASTA, if its ID exists in contig_dict,
    add contig_dict[contig_id]["bp"] = sequence_length
    """
    current_id = None
    current_len = 0

    def flush():
        nonlocal current_id, current_len
        if current_id is not None and current_id in contig_dict:
            contig_dict[current_id]["bp"] = current_len

    with open(fasta_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                flush()
                current_id = line[1:].split()[0]  # first token after ">"
                current_len = 0
            else:
                current_len += len(line)

        flush()

    return contig_dict


def main():
    BASE = Path(".")
    TXT_DIR = BASE / "cat_txt"
    FASTA_DIR = BASE / "fastas"

    all_rows = []

    for txt_file in TXT_DIR.glob("*.txt"):
        # infer expedition name
        expedition = txt_file.name.replace(".contig2classification.txt", "")
        fasta_file = FASTA_DIR / f"{expedition}.fasta"

        # build contig-level info
        contig_info = build_genus_dict_from_tsv(txt_file, threshold=THRESHOLD)
        contig_info = add_contig_lengths_from_fasta(fasta_file, contig_info)

        # contig_info -> DataFrame
        df = pd.DataFrame.from_dict(contig_info, orient="index")

        # sum bp per genus
        genus_bp = df.groupby("genus", as_index=False)["bp"].sum()

        # normalize to percentages
        total_bp = genus_bp["bp"].sum()
        genus_bp["percentage"] = genus_bp["bp"] / total_bp

        # limit to 3 significant figures
        genus_bp["percentage"] = genus_bp["percentage"].apply(
            lambda x: float(f"{x:.3g}")
        )

        # reshape to one row per expedition
        row = genus_bp.set_index("genus")["percentage"].to_dict()
        row["Expedition"] = expedition

        all_rows.append(row)

    # combine all expeditions
    final_df = pd.DataFrame(all_rows).fillna(0.0)

    # reorder columns: Expedition first, genera alphabetically
    cols = ["Expedition"] + sorted(c for c in final_df.columns if c != "Expedition")
    final_df = final_df[cols]

    # write output
    final_df.to_csv("bp_taxonomy.csv", sep=";", index=False)

if __name__ == "__main__":
    main()
