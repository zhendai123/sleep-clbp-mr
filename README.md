# sleep-clbp-mr

Analysis code for: **"Genetic liability to insomnia, but not obstructive sleep apnea, causally increases chronic low back pain: triangulating Mendelian randomization and NHANES evidence"**

## Contents
- `analysis_pipeline.py` - core MR pipeline: instrument harmonisation, IVW, MR-Egger,
  weighted median, Steiger directionality, two-step MR mediation functions.

## Data (all publicly available, no download included)
- NHANES 2009-2010: https://www.cdc.gov/nchs/nhanes (files DEMO_F, ARQ_F, SLQ_F, DPQ_F,
  MCQ_F, SMQ_F, ALQ_F, BMX_F, BPX_F, CRP_F, CBC_F)
- Insomnia GWAS: OpenGWAS ukb-b-3957
- CLBP / OSA / sciatica GWAS: FinnGen release 11 (M13_LOWBACKPAIN, G6_SLEEPAPNO,
  M13_LOWBACKPAIN_OR_AND_SCIATICA)
- Inflammatory proteins: GWAS Catalog GCST90274758-GCST90274848 (SCALLOP Olink panel)

## Usage
1. Download GWAS summary statistics from the sources above.
2. Extract instruments at genome-wide significance (P < 5e-8, F > 10), clump by linkage disequilibrium (r² < 0.001, 10 Mb window, 1000 Genomes European reference) via the OpenGWAS clumping API, and harmonise with outcome data as CSV files with the columns listed in analysis_pipeline.py. Distance-based clumping windows (500 kb, 1 Mb, 2 Mb) were evaluated as sensitivity analyses.
3. Run: `python analysis_pipeline.py`
