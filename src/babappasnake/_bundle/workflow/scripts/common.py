import math
import re
from pathlib import Path


STOP_CODONS = {"TAA", "TAG", "TGA"}


def ensure_parent(path):
    Path(path).parent.mkdir(parents=True, exist_ok=True)


def read_fasta_records(path):
    records = []
    header = None
    sequence_lines = []
    with open(path) as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append(
                        {
                            "id": header.split()[0],
                            "description": header,
                            "seq": "".join(sequence_lines),
                        }
                    )
                header = line[1:].strip()
                sequence_lines = []
            else:
                sequence_lines.append(line)
    if header is not None:
        records.append(
            {
                "id": header.split()[0],
                "description": header,
                "seq": "".join(sequence_lines),
            }
        )
    return records


def write_fasta_records(records, path, line_width=60):
    ensure_parent(path)
    with open(path, "w") as handle:
        for record in records:
            handle.write(f">{record['id']}\n")
            sequence = record["seq"]
            for index in range(0, len(sequence), line_width):
                handle.write(sequence[index:index + line_width] + "\n")


def sanitize_identifier(value):
    sanitized = re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("_")
    return sanitized or "item"


def canonicalize_taxon(value):
    return re.sub(r"\s+", " ", re.sub(r"[^A-Za-z0-9]+", " ", value)).strip().lower()


def guess_taxon_from_identifier(identifier):
    tokens = [token for token in re.split(r"[^A-Za-z]+", identifier) if token]
    for left, right in zip(tokens, tokens[1:]):
        if re.fullmatch(r"[A-Z][a-z]+", left) and re.fullmatch(r"[a-z][a-z0-9-]*", right):
            return f"{left} {right}"
    return None


def has_internal_stop_codon(sequence):
    ungapped = sequence.upper().replace("U", "T").replace("-", "")
    if len(ungapped) < 6:
        return False
    for index in range(0, len(ungapped) - 3, 3):
        codon = ungapped[index:index + 3]
        if len(codon) < 3 or "N" in codon or "?" in codon:
            continue
        if codon in STOP_CODONS:
            return True
    return False


def benjamini_hochberg(p_values):
    indexed = sorted(enumerate(p_values), key=lambda item: item[1])
    corrected = [1.0] * len(p_values)
    previous = 1.0
    total = len(p_values)
    for rank, (original_index, p_value) in enumerate(reversed(indexed), start=1):
        adjusted = min(previous, p_value * total / (total - rank + 1))
        corrected[original_index] = adjusted
        previous = adjusted
    return corrected


def mixed_chi_square_pvalue(lrt_statistic):
    if lrt_statistic <= 0:
        return 1.0
    # 50:50 mixture of a point mass at 0 and chi-square(df=1).
    from scipy.stats import chi2

    return 0.5 * chi2.sf(lrt_statistic, df=1)


def make_unique_member_ids(candidate_labels):
    seen = {}
    resolved = []
    for label in candidate_labels:
        base = sanitize_identifier(label)
        count = seen.get(base, 0)
        if count == 0:
            resolved.append(base)
        else:
            resolved.append(f"{base}__{count + 1}")
        seen[base] = count + 1
    return resolved


def relative_path(path, start):
    try:
        return str(Path(path).resolve().relative_to(Path(start).resolve()))
    except ValueError:
        return str(path)


def as_bool(value):
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "yes", "y"}


def retention_cutoff(max_count, min_ratio, max_loss):
    return max(int(math.ceil(max_count * min_ratio)), int(max_count - max_loss))
