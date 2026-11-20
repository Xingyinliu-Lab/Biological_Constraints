import csv
import hashlib
import json
import math
import pathlib
import shutil
import subprocess
import sys
import tempfile

SCORING_DIR = pathlib.Path(__file__).resolve().parents[2] / "scoring"
sys.path.insert(0, str(SCORING_DIR))
from ps10_refset import validate_records


RECORD_SHA = "c988e89f50d92ecdda8805ba82ae7f9f95483eea5bd342bdcedfd7e41c9a0dff"
COUNT = 2031
CONFIG = {"ct": "N", "nt": "N", "ph": 7.0, "te_kelvin": 298, "io_molar": 0.1, "scalar": "Aggregation"}


def sha(path):
    h = hashlib.sha256()
    with pathlib.Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main(records_path, binary, output_dir):
    records_path, binary, output_dir = map(pathlib.Path, (records_path, binary, output_dir))
    actual_records_sha256 = sha(records_path)
    rows = validate_records(records_path)
    ids = [row["record_id"] for row in rows]
    if len(rows) != COUNT or len(set(ids)) != COUNT or any(len(x) >= 25 for x in ids):
        raise RuntimeError("record population or TANGO name limit mismatch")
    if output_dir.exists():
        raise RuntimeError("refusing to overwrite existing TANGO output")
    output_dir.mkdir(parents=True)
    scores = {}
    for batch_no, start in enumerate(range(0, len(rows), 900), 1):
        batch = rows[start : start + 900]
        with tempfile.TemporaryDirectory(prefix="ps10_tango_") as td:
            td = pathlib.Path(td)
            name = "input" + str(batch_no) + ".in"
            path = td / name
            with path.open("w", encoding="ascii", newline="\n") as handle:
                for row in batch:
                    handle.write(f"{row['record_id']} {CONFIG['ct']} {CONFIG['nt']} {CONFIG['ph']} {CONFIG['te_kelvin']} {CONFIG['io_molar']} {row['normalized_peptide_identity']}\n")
            subprocess.run([str(binary)], cwd=td, input=f"N\n{name}\n", text=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, check=True)
            aggregates = list(td.glob("*_aggregation.*"))
            if len(aggregates) != 1:
                raise RuntimeError("missing TANGO aggregate output")
            aggregate = aggregates[0]
            parsed = list(csv.DictReader(aggregate.open(encoding="utf-8"), delimiter="\t"))
            if len(parsed) != len(batch):
                raise RuntimeError("TANGO batch output count mismatch")
            for item in parsed:
                key = item["Sequence"].strip()
                value = float(item[CONFIG["scalar"]])
                if key in scores or not math.isfinite(value):
                    raise RuntimeError("invalid or duplicate TANGO score")
                scores[key] = value
    if set(scores) != set(ids):
        raise RuntimeError("TANGO score population closure failure")
    score_path = output_dir / "TANGO_SCORES.tsv"
    with score_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["record_id", "tango_aggregation"], delimiter="\t")
        writer.writeheader()
        writer.writerows({"record_id": record_id, "tango_aggregation": repr(scores[record_id])} for record_id in ids)
    report = {"status": "PASS", "record_count": len(ids), "records_sha256": actual_records_sha256, "registered_records_sha256": RECORD_SHA, "registered_records_bytes_match": actual_records_sha256 == RECORD_SHA, "binary_sha256": sha(binary), "scores_sha256": sha(score_path), "configuration": CONFIG, "batch_count": math.ceil(len(rows) / 900)}
    (output_dir / "TANGO_FORWARD_REPORT.json").write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main(*sys.argv[1:])
