# Predicting HER2/EGFR selectivity from ligand structure

**MSBI-565 Computational Biomedicine and Drug Discovery — Final Project**
Instructor: Abdullah Al Mamun

A QSAR study asking whether ligand structure predicts not just how potently a
compound inhibits HER2, but **how selective** it is for HER2 over its closest
relative EGFR — and whether measurement noise, rather than model choice, is
what limits the answer.

---

## 1. The question

HER2 (ERBB2) and EGFR (ERBB1) have near-identical ATP-binding sites, and
clinical ErbB inhibitors span the whole selectivity range: lapatinib is dual,
neratinib is a covalent pan-HER inhibitor, tucatinib is deliberately
HER2-selective to avoid the rash and diarrhoea caused by EGFR inhibition in
skin and gut. Selectivity is therefore a design objective, not a footnote.

Three models, all evaluated by identical code so the comparison is fair:

| | Target | Trained on |
|---|---|---|
| **Model A** | `pHER2` | all curated HER2 compounds |
| **Model B** | `pEGFR` | all curated EGFR compounds |
| **Model C** | `Δ = pHER2 − pEGFR` | compounds measured against both |

Model C is compared against the naive **baseline of subtracting Model A's and
Model B's predictions**, and both are compared against a ceiling derived from
the reproducibility of the underlying measurements.

---

## 2. The central result, stated up front

Compounds in this dataset that were assayed for HER2 in **more than one
independent study** disagree with themselves by a pooled standard deviation of
**0.63 log units**, and 30% of them disagree by more than a full log unit.
That is the noise in the labels, and it caps any model at

```
R²_max  =  1 − s_rep² / s_total²  ≈  0.69
RMSE_min ≈ 0.63 log units  (a factor of ~4 in IC50)
```

Selectivity is a **difference** of two such measurements, so its noise is
`√2 × 0.63 ≈ 0.89` log units — while its dynamic range is narrower than that
of potency, because most ErbB inhibitors hit both targets. Noise up, signal
down. The project's quantitative prediction is that the selectivity ceiling
sits far below 0.69, and that neither Model C nor the baseline can exceed it.

Any QSAR paper reporting R² above 0.69 on this kind of data is almost
certainly measuring analog leakage, not generalisation.

---

## 3. Install and run

```bash
git clone <your-repo-url>
cd her2-egfr-qsar
pip install -r requirements.txt
```

Run from the **repository root**, in order:

```bash
python src/01_explore_her2.py       # offline. Curation + noise floor
python src/02_fetch_egfr.py         # NEEDS INTERNET. Downloads EGFR
python src/03_curate.py             # NEEDS RDKit. Standardise + matched set
python src/04_features.py           # NEEDS RDKit. Descriptors + scaffolds
python src/05_models.py             # Models A and B + sanity checks
python src/06_selectivity.py        # Model C, baseline, AD, drivers
```

Step 02 pauses after building the assay catalogue:

```bash
python src/02_fetch_egfr.py --catalogue-only   # writes data/raw/egfr_assay_catalogue.csv
#   ... open it, review the `keep` and `drop_reason` columns, edit as needed ...
python src/02_fetch_egfr.py --download-only    # then download the kept assays
```

**Review that catalogue.** The auto-triage is a first pass. It drops mutant
EGFR assays (`L858R`, `T790M`, `E746-A750del` …), cell-proliferation and
viability assays, and selectivity/ratio endpoints — each for a stated reason
in `drop_reason`. Pooling a kinase-domain mutant into the wild-type EGFR
endpoint would make Δ meaningless, since gefitinib is ~100× more potent on
L858R than on wild type.

### Running in Google Colab

RDKit and network access both work out of the box:

```python
!pip install rdkit
!git clone <your-repo-url> && cd her2-egfr-qsar && python src/01_explore_her2.py
```

### Reproducibility

Every split, model and shuffle derives from a single seed, `SEED = 20261109`,
set in `src/common.py`. Re-running any script reproduces the tables exactly.

---

## 4. Layout

```
her2-egfr-qsar/
├── README.md
├── requirements.txt
├── data/
│   ├── raw/                     # downloaded, never edited by code
│   │   ├── combined_pubchem_data.csv            (HER2, supplied)
│   │   ├── combined_pubchem_data_assay_catalogue.csv
│   │   ├── egfr_assay_catalogue.csv             (from 02)
│   │   └── egfr_pubchem_data.csv                (from 02)
│   └── processed/               # everything generated
├── src/
│   ├── common.py                # paths, seed, constants, plot style
│   ├── modeling.py              # splits, model zoo, evaluation loop
│   ├── 01_explore_her2.py       # curation waterfall + noise floor
│   ├── 02_fetch_egfr.py         # PubChem download + assay triage
│   ├── 03_curate.py             # standardisation + matched set
│   ├── 04_features.py           # 3 feature sets + Murcko scaffolds
│   ├── 05_models.py             # Models A and B
│   └── 06_selectivity.py        # Model C, baseline, AD, drivers
├── notebooks/
│   └── walkthrough.ipynb        # narrated end-to-end run
├── results/
│   ├── figures/                 # fig1 … fig10, referenced by the report
│   └── tables/                  # table1 … table10
└── report/
    └── report.md                # the 5-page write-up
```

### What each script produces

| Script | Tables | Figures |
|---|---|---|
| `01_explore_her2.py` | 1 waterfall, 2 endpoints, 3 **noise floor** | 1–4 |
| `03_curate.py` | 4 matched-set summary | 5 selectivity landscape |
| `04_features.py` | — | 6 chemical space |
| `05_models.py` | 5 Model A, 6 Model B, 7 sanity checks | 7–8 |
| `06_selectivity.py` | 8 Model C, 9 **C vs baseline**, 10 drivers | 9–10 |

---

## 5. Method decisions, and why

**Endpoint.** `pActivity = −log₁₀(molar)`, pooling IC50/Ki/Kd/EC50 in nM. The
original type is retained per row so the pooling can be revisited.

**Regression, not classification.** The HER2 set has 1,294 `Active`, 900
`Unspecified`, 52 `Inconclusive` and **4 `Inactive`** outcomes. Four negatives
cannot support an honest classifier. The potency value is the real signal, so
we regress on it.

**Censored values excluded.** 337 rows report `>` and 21 report `<`. A value
of "> 10 µM" means inactive up to the highest concentration tested — it is not
10 µM. Keeping them as exact values would compress the inactive tail and
flatter every model. They are quantified in Figure 2b and set aside.

**Standardisation.** Largest organic fragment → neutralise charges → canonical
isomeric SMILES. Identical for both targets, because any asymmetry would break
the join and silently shrink the matched set. Without neutralisation a free
base and its hydrochloride become different compounds, the replicate is
missed, and the pair can land on opposite sides of a split.

**Scaffold splitting, always.** `GroupShuffleSplit` on Bemis-Murcko scaffolds,
5 splits, 25% of *scaffolds* held out. Random splits are reported too, but only
as the deliberately optimistic contrast; the gap between them is the size of
the analog-leakage effect.

**Uncertainty on everything.** Every metric is mean ± SD across the 5 splits.
A single split of ~1,000 compounds has a standard error large enough to
reverse model rankings.

**Leakage control in the baseline.** For each selectivity split, Models A and
B are retrained from scratch with **every scaffold in the test fold removed**
from their training sets. Skipping this is the easiest way to produce a
beautiful, meaningless baseline.

**Sanity checks.** A `DummyRegressor` in every table; y-scrambling over 20
rounds (R² must collapse to ≈ 0); and the random-vs-scaffold contrast.

**No hyperparameter tuning.** Models use scikit-learn defaults with minor
documented exceptions. Tuning inside the evaluation loop without a nested CV
would leak the test fold, and a proper nested CV over 6 models × 3 feature
sets is beyond the project's compute budget. This is stated as a limitation
rather than hidden.

---

## 6. If the matched set is too small

The one real risk: Model C needs compounds measured against **both** targets.
HER2 and EGFR are co-assayed constantly in the same medicinal-chemistry
papers, so the overlap should be substantial — but `03_curate.py` prints the
number and warns below `MIN_COMPOUNDS = 150`.

If you land under ~150, two fallbacks, in order of preference:

1. **Relax the join.** Allow censored EGFR values as a *lower bound* on Δ and
   model selectivity as a classification (HER2-selective / dual /
   EGFR-selective) rather than a regression. Fewer assumptions, more data.
2. **Drop Model C, keep the comparison.** Report Models A and B plus the
   noise-floor analysis, and present the A−B baseline as the selectivity
   estimate with its propagated error bar. The noise argument stands on its
   own and is still the most interesting part of the project.

---

## 7. Data sources

All bioactivity data is from **PubChem BioAssay**, which aggregates
depositions from ChEMBL and from individual publications.

- HER2 / ERBB2 — NCBI Gene ID **2064**, 253 curated assays (supplied with the
  course materials as `combined_pubchem_data.csv`)
- EGFR / ERBB1 — NCBI Gene ID **1956**, downloaded by `src/02_fetch_egfr.py`

PubChem PUG REST: `https://pubchem.ncbi.nlm.nih.gov/rest/pug`. Requests are
rate-limited to 4/second, with exponential backoff on 503.

Raw downloads are **never** modified in place; every transformation writes to
`data/processed/`.

---

## 8. Known limitations

1. **Pooled assay conditions.** Measurements come from ~138 HER2 assays run by
   different groups at different ATP concentrations. ATP-competitive inhibitor
   IC50 depends on [ATP], and we do not correct for it. This is a large part
   of the 0.63 log-unit noise floor.
2. **Pooled endpoint types.** IC50, Ki, Kd and EC50 are not the same quantity.
   IC50 dominates (1,385 of 1,440 measurements), so the effect is small, but
   it is real.
3. **No hyperparameter tuning** (see above).
4. **2D descriptors only.** No conformers, no docking, no explicit modelling
   of the ATP pocket — so the model cannot reason about *why* a substituent is
   tolerated by one kinase and not the other.
5. **Publication bias.** PubChem over-represents potent compounds; weak and
   inactive molecules are under-reported, which truncates the lower tail.
6. **Median aggregation** of replicate measurements discards the spread. A
   weighted fit using per-compound SD would use the data better.
