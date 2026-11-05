# Exact upper bounds for three binary constant-weight codes

A(*n,d,w*) is the largest size of a binary code of length *n*, minimum Hamming distance *d* and constant weight *w*. Each rational bound below is an upper bound on the size of every such code. Taking its floor gives the integer upper bound.

| (n,d,w) | Previous integer bound | Exact rational upper bound | Exact decimal | Integer upper bound |
| --- | --- | --- | --- | --- |
| (18,6,7) | 312 | 3868791767977833258518625959723/12500000000000000000000000000 | 309.50334143822666068149007677784 | 309 |
| (19,6,7) | 463 | 46047385392564732804912812389921/100000000000000000000000000000 | 460.47385392564732804912812389921 | 460 |
| (19,6,8) | 673 | 16683417296961520597884426836631/25000000000000000000000000000 | 667.33669187846082391537707346524 | 667 |

The previous upper bounds appear in the distance-six section of [Brouwer’s table](https://aeb.win.tue.nl/codes/Andw.html).

## Paper

The [paper (PDF)](https://github.com/wille614/cw-code-upper-bounds/releases/download/v1.2.0/paper.pdf) gives all constructions and proofs; the [LaTeX source](https://github.com/wille614/cw-code-upper-bounds/blob/v1.2.0/paper/paper.tex) is also available.

## Certificates and models

Each augmented four-point SDP uses anticode vacancy, occupancy and cardinality inequalities. The cardinality premise is obtained by flooring an exact Polak four-point bound. The three baseline certificates establish the premises 312, 463 and 673, respectively. The package contains the baseline and augmented SDPs and their rational dual certificates. The derivation of every constraint family, and the proof that every code gives a feasible point of the augmented programs, is in [proofs/PROOFS.md](proofs/PROOFS.md).

The release models have normalized metadata and exactly the same mathematical content as the original computations, and the certificates have been rebound to them and verified.

Each certificate records rational dual matrices, their LDLᵀ factors, residual compensation and hashes of its JSON and SDPA models. A model’s `cardinality_provenance.artifact` is relative to the top of this repository. Its adjacent SHA-256 identifies the exact baseline certificate. The baseline certificate for A(19,6,8) ≤ 673 is byte-identical to the [published certificate](https://github.com/wille614/a4-19-6-8-certificate/releases/tag/v1.0.0).

## Reproduce the bounds offline

1. Start in the directory containing this README. Use CPython 3.11 or later with the standard library and assertions enabled.
2. Download the 18 gzip data assets from [Release v1.0.0](https://github.com/wille614/cw-code-upper-bounds/releases/tag/v1.0.0) and place them in the directories indicated below. The data and certificates are unchanged in v1.2.0. In an assembled copy they are already in place. Each filename is listed in its directory’s checksum list.

   | Directory | Files per point |
   | --- | --- |
   | [baselines/](baselines/SHA256SUMS) | `B<n>-<d>-<w>.certificate.json.gz`, `.model.json.gz` and `.dat-s.gz` |
   | [certificates/](certificates/SHA256SUMS) | `A<n>-<d>-<w>.certificate.json.gz` |
   | [data/](data/SHA256SUMS) | `A<n>-<d>-<w>.model.json.gz` and `.dat-s.gz` |

3. Run the commands for your platform. Choose a fresh replay directory. Preparation verifies the six asset hashes per point, decompresses both models for each SDP, and checks the augmented model’s baseline reference, rational bound and floor. The two verifier invocations independently establish the premise and the improved bound.

**Windows PowerShell**

```powershell
$ErrorActionPreference = 'Stop'
foreach ($point in '18-6-7','19-6-7','19-6-8') {
    py -3 -I prepare_inputs.py $point "../cw-replay/$point"
    if ($LASTEXITCODE -ne 0) { throw "Input preparation failed for $point" }
    foreach ($group in 'baseline','upper') {
        py -3 -I verifier/verify_certificate.py "../cw-replay/$point/$group"
        if ($LASTEXITCODE -ne 0) { throw "Verification failed for $point ($group)" }
    }
}
```

**Linux or macOS**

```sh
set -eu
for point in 18-6-7 19-6-7 19-6-8
do
  python3 -I prepare_inputs.py "$point" "../cw-replay/$point"
  for group in baseline upper
  do
    python3 -I verifier/verify_certificate.py "../cw-replay/$point/$group"
  done
done
```

Each verifier prints a JSON object and writes `verification.json` in its input directory. Expect six `PASS` results with exit code 0. The integer-bound pairs (baseline, augmented) are (312, 309), (463, 460) and (673, 667). The three augmented runs check 46677, 52716 and 116007 blocks, covering 13834807, 15981482 and 66812277 SDPA coefficients, respectively. Check one point at a time.

## Attribution and licenses

Copyright © 2026 Jiang Chang. Code is licensed under the [Apache License 2.0](LICENSE), and certificates and model data under [CC BY 4.0](LICENSE-DATA). Cite this package using [CITATION.cff](CITATION.cff). The unmodified exact verifier comes from [Jiang Chang’s four-point certificate repository](https://github.com/wille614/a4-19-6-8-certificate/tree/79dbdd79e4c9a4e24e54375d2c4930cab45bd708).
