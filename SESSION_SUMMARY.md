# AlphaGenome QTL Validation - My Development Notes

**Author**: George Stephenson  
**Lab**: LAYER Lab, CU Boulder

## My Project Goal

I'm validating AlphaGenome's regulatory variant effect predictions using chromatin accessibility QTLs (caQTLs) and histone modification QTLs (hQTLs).

## Critical Issues Resolved

### Issue #1: Placeholder Alleles (Week of Debugging!)

My previous approach (`alphagenome_phase3_qtl_benchmark`) used placeholder A→G alleles for all 15,647 variants, resulting in r≈0.03 after a WEEK of failed attempts.

**What I did wrong**: Scored wrong molecular changes (e.g., T→G instead of actual T→C variant)

**How I fixed it**: Fetch real alleles from dbSNP using rs IDs in QTL summary statistics

**Result**: r=0.40 for caQTLs (13x improvement!)

### Issue #2: Directional Effects

Even with real alleles, my initial tests still showed r≈0 due to using absolute value correlation.

**What I did wrong**: Used `|quantile_score|` vs `beta`, destroying directional signal

**How I fixed it**: Use raw quantile_score vs beta (both are directional effects)

**Result**: Proper correlation capturing increase/decrease in chromatin accessibility

## Validated Datasets

### caQTLs (Chromatin Accessibility)
- Source: Nedelec et al., GSE86886, CD4+ T cells
- Variants: 3,271 with real alleles (98.7% success rate)
- Modalities: ATAC-seq, DNase-seq, H3K27ac ChIP-seq
- Tissue filter: CD4+ T cells (CL:0000624)
- **Results**: r=0.4032 (Spearman), r=0.4089 (Pearson), p < 10⁻¹²⁸

### hQTLs (Histone Modifications)  
- Source: Pelikan et al., GSE116193, lymphoblastoid cell lines
- Variants: 6,259 with real alleles (99.97% success rate)
- Modalities: H3K27ac and/or H3K4me1 ChIP-seq
- Tissue filter: B cells (CL:0000236)
- **Results**: In progress... ⏳

### dsQTLs (DNA Shape) - FAILED
- Source: GSE31388 (2012 dataset)
- Attempted: 6,070 variants
- **Result**: 0% success - hg19 positions not in modern dbSNP
- **Lesson**: Old datasets (pre-2015) often lack annotations needed for allele fetching
- **Status**: Skipped, but script still supports it if I find better data

#### Step 2: Score with AlphaGenome (`02_predict.py`)
- Follows official `batch_variant_scoring.ipynb` pattern exactly
## Implementation Details

### Allele Fetching
- API: myvariant.info with hg38 assembly
- Batch size: 100 variants per request
- Rate limiting: 0.5 sec between batches
- Handles multi-allelic sites (takes first alt allele)

### Tissue-Matched Scoring
- caQTLs: CD4+ T cells (CL:0000624) - matches Nedelec et al. source
- hQTLs: B cells (CL:0000236) - LCLs are immortalized B cells
- Filters AlphaGenome predictions to tissue-specific tracks only

### Modality Handling
- Deduplicates scorers for multi-modality variants
- Example: H3K27ac,H3K4me1 → single CHIP_HISTONE scorer (not duplicate)
- Groups variants by modality before scoring

### Evaluation
- Directional correlation: quantile_score vs beta (both signed effects)
- Reports: Spearman r (rank), Pearson r (linear), p-values

## Key Lessons

1. **Always verify input data matches model expectations**
   - I wasted a WEEK using placeholder alleles that produced meaningless results
   - Always fetch real alleles from reference databases

2. **Understand your evaluation metrics**
   - QTL betas are directional (increase/decrease chromatin state)
   - I mistakenly used absolute value which destroyed biological signal

3. **Tissue matching is critical**
   - Regulatory effects are highly tissue-specific
   - I must filter to source cell type for accurate predictions

4. **Handle multi-modality variants correctly**
   - Some variants affect multiple related marks (H3K27ac + H3K4me1)
   - I had to deduplicate scorers to avoid API errors

5. **Systematic debugging with increasing sample sizes**
   - Small samples (n=5) can show spurious correlations
   - I now test on n=20, 100 before running full dataset

6. **Read official examples carefully**
   - batch_variant_scoring.ipynb showed me the correct pattern
   - Tissue filtering turned out to be critical

7. **Simplify when stuck**
   - My ISM approach was too complex
   - Real alleles from API was the simpler solution

8. **Old datasets can be problematic**
   - 2012 dsQTL dataset lacked modern dbSNP annotations
   - Pre-2015 data often can't be used with current pipelines

## My Next Steps

1. ✅ caQTLs complete: r=0.40 (SUCCESS!)
2. ⏳ hQTLs running in background (~1.5 hours)
3. ❌ dsQTLs failed (old dataset, skipped)
4. Generate final plots and benchmark report

## Directory Structure

```
alphagenome-qtl-validation/
├── data/
│   ├── raw/               # Original QTL files  
│   └── processed/         # With REAL alleles from API
├── scripts/
│   ├── 01_prepare_qtls.py    # Fetch real alleles
│   ├── 02_predict.py          # Score with AlphaGenome
│   └── 03_evaluate.py         # Compute correlations
├── results/
│   ├── predictions/       # AlphaGenome scores
│   └── tables/           # Evaluation metrics
└── logs/                 # Execution logs
```

## My Commands

```bash
# Prepare data with real alleles
python scripts/01_prepare_qtls.py --datasets caQTLs hQTLs

# Run predictions
python scripts/02_predict.py --datasets caQTLs hQTLs --force

# Evaluate  
python scripts/03_evaluate.py --datasets caQTLs
```

---
**Date**: November 17, 2025  
**Status**: Full run in progress, test results promising (r=0.64)
