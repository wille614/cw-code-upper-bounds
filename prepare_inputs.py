"""Prepare both SDPs for one bound and check the reference to its four-point certificate."""
import argparse
from fractions import Fraction
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import tempfile


ROOT = Path(__file__).resolve().parent
POINTS = ("18-6-7", "19-6-7", "19-6-8")


def check_asset(path):
    manifest = (path.parent / "SHA256SUMS").read_text(encoding="utf-8")
    expected = {name: digest for digest, name in (line.split() for line in manifest.splitlines())}
    with path.open("rb") as stream:
        actual = hashlib.file_digest(stream, "sha256").hexdigest()
    if actual != expected[path.name]:
        raise ValueError(f"SHA-256 mismatch for {path.name}")
    return actual


def prepare(folder, cert, model, sdpa):
    folder.mkdir()
    cert_hash = check_asset(cert)
    shutil.copyfile(cert, folder / "certificate.json.gz")
    for archive, filename in ((model, "compressed_model.json"), (sdpa, "compressed.dat-s")):
        check_asset(archive)
        with gzip.open(archive, "rb") as src, (folder / filename).open("xb") as dst:
            shutil.copyfileobj(src, dst, 1024 * 1024)
    return cert_hash


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("point", choices=POINTS)
    parser.add_argument("output", type=Path, help="new directory for the two verifier inputs")
    args = parser.parse_args()
    output = args.output.resolve()
    if output.is_relative_to(ROOT) or output.exists():
        parser.error("choose a new output directory outside this checkout")
    point = args.point
    base = ROOT / "baselines"
    cert = base / f"B{point}.certificate.json.gz"

    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=output.parent) as temp:
        work = Path(temp) / "inputs"
        work.mkdir()
        base_hash = prepare(work / "baseline", cert, base / f"B{point}.model.json.gz",
                            base / f"B{point}.dat-s.gz")
        prepare(work / "upper", ROOT / "certificates" / f"A{point}.certificate.json.gz",
                ROOT / "data" / f"A{point}.model.json.gz", ROOT / "data" / f"A{point}.dat-s.gz")
        with gzip.open(work / "baseline" / "certificate.json.gz", "rt", encoding="utf-8") as stream:
            bound = Fraction(json.load(stream)["bound"])
        with (work / "upper" / "compressed_model.json").open(encoding="utf-8") as stream:
            ref = json.load(stream)["cardinality_provenance"]
        if (ref["artifact"] != cert.relative_to(ROOT).as_posix()
                or ref["certificate_sha256"] != base_hash
                or Fraction(ref["rational_upper"]) != bound
                or ref["bound"] != bound.numerator // bound.denominator):
            raise ValueError("The model's cardinality cap does not match its four-point certificate")
        work.rename(output)
    print(f"Prepared {output} (four-point cap {ref['bound']})")


if __name__ == "__main__":
    main()
