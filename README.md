# Type-D infinitesimal braid deformations: reproducibility materials

Companion code and exact certificates for **A Uniform Calculation of Weight-Preserving Infinitesimal Braid Deformations in Type D**, contributed by **Hangang Wang** with extensive OpenAI ChatGPT and Codex assistance.

## Scope and status

The accompanying manuscript studies the full classical-weight-preserving linearized braid equation at the standard type-$D_n$ vector braiding over $K=\mathbb{Q}(q)$, with $q$ indeterminate. Its uniform theorem states

$$\dim_K \ker(\mathcal{L}_Q|_{\mathcal{W}_n})=\frac{(n+1)(n+2)}2,\qquad n\ge4.$$

The written proof supplies the all-rank statement. These programs check transcription, implementations, finite-rank minors and the known-family parameter dictionary. More finite-rank checks do not constitute a proof for every $n$. The explicit integrable family is already part of earlier multiparameter constructions. Bibliographic priority for the full tangent-space exhaustion statement remains unresolved.

Multiple AI-assisted audits have found no counterexample. No independent human-specialist verification or proof-assistant formalization is claimed. No conclusion at every special value of $q$, unrestricted deformation classification, or global classification is claimed.

## Run

Python **3.10 or later** is required. From this repository's root:

```sh
python code/run_all.py
```

The default checks use only Python's standard library. For the separate symbolic audit:

```sh
python -m pip install -r requirements-audit.txt
python code/run_all.py --include-independent-audit
```

Results are written into `results/`; reruns update those local JSON files. Timing fields can vary. Compare mathematical values, ranks, hashes, counts and PASS statuses rather than expecting byte-identical timing fields. Native MATLAB is not used by this suite and was not rerun for the October verification.

## Files and verification boundaries

| Directory | Contents |
|---|---|
| `code/` | Current portable programs and exact Laurent arithmetic |
| `independent_review/` | Separately implemented AI-assisted coefficient checks |
| `input_certificates/` | Six archived modular pivot certificates for $D_4$–$D_7$ |
| `results/` | Recorded outputs of the replayed checks |
| `AI_AND_PROVENANCE.md` | Historical contributions, AI use, and limits of verification |
| `MANIFEST.json` | SHA-256 inventory of this snapshot, excluding the inventory itself |

The independent replay comprises **150 + 628 modular comparisons**, **28 additional row-$N$ comparisons**, and **5 exact symbolic equality assertions**. The main finite-rank certificates use $q=2$ and $p=1000000007$; alternate $D_4,D_5$ certificates use $q=3$, $p=1000000009$. The dictionary check compares the complete operator coefficients and twist constraints for $n=4,5$.

Historical comments may retain equation numbers from earlier drafts. Identify an equation by its formula and tensor input/output. In the formatted manuscript the row-$N$ equation for input $kN1$ and output $Nk1$ is (4.7).

## Manuscript and citation

The manuscript's public Hexagon identifier will be added after publication. No placeholder identifier or DOI is asserted here. `CITATION.cff` describes this code snapshot. Cite the **specific commit used**, rather than only the moving `main` branch; the manuscript records that commit once this snapshot is uploaded.

January 2026 source dates refer to the historical computations, not the later all-rank theorem. This curated public release is dated **27 October 2026**. Private correspondence, third-party textbooks and complete private handoff archives are excluded.

## Related literature

The code is a companion to a manuscript that compares the family with [Isaev–Ogievetsky–Pyatov](https://arxiv.org/abs/math/0509251v1) and [Martin–Tsymbaliuk](https://arxiv.org/abs/2508.10882v1), among other sources. These links are references, not redistributed third-party source material. They do not certify novelty of the exhaustion theorem.

## License

No repository license has yet been selected. Public availability is not a declaration that the code has an open-source license. The manuscript's eventual Hexagon license is selected separately by the contributor.
