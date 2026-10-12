# Sleep disturbance and chronic low back pain - analysis code

Companion code for: "Sleep disturbance and chronic low back pain: triangulated evidence
from a nationally representative survey and Mendelian randomization"

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
2. Extract instruments (P < 5e-8, clumping window 500 kb, F > 10) and harmonise with
   outcome data as CSV files with the columns listed in `analysis_pipeline.py`.
3. Run: `python analysis_pipeline.py`
