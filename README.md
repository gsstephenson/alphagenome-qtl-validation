# AlphaGenome QTL Validation

**Author**: George Stephenson  
**Lab**: LAYER Lab, CU Boulder

I validate AlphaGenome's regulatory variant effect predictions using chromatin accessibility QTLs (caQTLs) and histone modification QTLs (hQTLs).

## Overview

I validate AlphaGenome by correlating its predicted variant effects with experimentally measured QTL effect sizes. This tests whether AlphaGenome can identify regulatory variants that affect chromatin state.

**My Approach**: 
1. **Fetch real alleles from dbSNP** using rs IDs from QTL summary statistics
2. **Score actual variants** with AlphaGenome using correct ref/alt bases
3. **Tissue-matched predictions** - filter to source cell types (CD4+ T cells, B cells)
4. **Directional correlation** - compare quantile_score vs beta (both directional effects)

## Results

### caQTLs (Chromatin Accessibility)
- **Dataset**: 3,271 variants from CD4+ T cells (Nedelec et al., GSE86886)
- **Spearman r = 0.4032** (p < 10⁻¹²⁸) ✅
- **Pearson r = 0.4089** (p < 10⁻¹²⁸) ✅
- **Tissue**: CD4+ T cells (CL:0000624)

### hQTLs (Histone Modifications)
- **Dataset**: 6,259 variants from B cells/LCLs (Pelikan et al., GSE116193)
- **Status**: Running predictions... ⏳
- **Modalities**: H3K27ac and H3K4me1 (CHIP_HISTONE scorer)
- **Tissue**: B cells (CL:0000236)

### dsQTLs (DNA Shape)
- **Status**: Skipped - 2012 dataset (GSE31388) lacks dbSNP annotations ❌
- **Issue**: hg19 positions not in modern databases, can't fetch real alleles
- **Note**: Script supports dsQTLs if better annotated data becomes available

These moderate-strong correlations demonstrate that AlphaGenome captures real regulatory variant effects!

## Pipeline

### Step 1: Fetch Real Alleles
```bash
python scripts/01_prepare_qtls.py --datasets caQTLs hQTLs
```
- Extracts rs IDs from QTL summary statistics
- Fetches real ref/alt alleles from myvariant.info API (batch size 100)
- Success rates:
  - **caQTLs**: 98.7% (3,274/3,317)
  - **hQTLs**: 99.97% (6,259/6,261)
  - **dsQTLs**: 0% (old dataset, no dbSNP annotations - skipped)
- Output: `data/processed/{dataset}.parquet` with real alleles

### Step 2: Score Variants with AlphaGenome
```bash
python scripts/02_predict.py --datasets caQTLs hQTLs --force
```
- Uses RECOMMENDED_VARIANT_SCORERS for each modality:
  - caQTLs: ATAC, DNASE, CHIP_HISTONE
  - hQTLs: CHIP_HISTONE (H3K27ac, H3K4me1)
- Filters to tissue-specific tracks matching source:
  - caQTLs: CD4+ T cells (CL:0000624)
  - hQTLs: B cells (CL:0000236)
- Aggregates across tracks: mean quantile_score per variant
- Speed: ~1.4 variants/sec
- Output: `results/predictions/{dataset}.parquet`

### Step 3: Evaluate Correlations
```bash
python scripts/03_evaluate.py --datasets caQTLs hQTLs
```
- Correlates **quantile_score vs beta** (directional effects)
- Computes: Spearman r, Pearson r, p-values
- Output: `results/tables/{dataset}_metrics.txt`

### Step 4: Generate Plots
```bash
python scripts/04_plots.py --datasets caQTLs hQTLs
```
- Correlation scatter plot
- Score and beta distributions
- Residual analysis
- Top/bottom variant effects
- Output: `results/plots/{dataset}_*.png`

## Key Implementation Details

**Allele Fetching:**
- API: `myvariant.info` with hg38 assembly
- Fields: `dbsnp.ref`, `dbsnp.alt`
- Rate limiting: 0.5 sec between batches
- Handles multi-allelic sites (takes first alt allele)

**Tissue-Matched Predictions:**
- caQTLs: Filter to CD4+ T cells (CL:0000624) - matches source tissue
- hQTLs: Filter to B cells (CL:0000236) - LCLs are immortalized B cells
- Critical for accurate tissue-specific regulatory predictions!

**Modality Handling:**
- Deduplicates scorers when variants have multiple modalities
- Example: H3K27ac,H3K4me1 → single CHIP_HISTONE scorer (not duplicate)
- Groups variants by modality and scores each group appropriately

**Evaluation Metric:**
- Uses RAW quantile_score (not absolute value!)
- Both scores and betas are directional effects
- Directional correlation captures sign of effect (increase/decrease accessibility)

## Data Sources

- **caQTLs**: Chromatin accessibility QTLs from CD4+ T cells
  - Source: Nedelec et al., GSE86886
  - Modalities: ATAC-seq, DNase-seq, H3K27ac ChIP-seq
  - Sample: 3,271 variants with effect sizes
  
- **hQTLs**: Histone modification QTLs from lymphoblastoid cell lines (LCLs)
  - Source: Pelikan et al., GSE116193
  - Modalities: H3K27ac and/or H3K4me1 ChIP-seq
  - Sample: 6,259 variants (3,502 both marks, 2,347 H3K27ac only, 410 H3K4me1 only)

- **dsQTLs**: DNA shape QTLs (SKIPPED)
  - Source: GSE31388 (2012 dataset)
  - Issue: hg19 positions not in modern dbSNP, can't fetch real alleles
  - Note: Script still includes dsQTL loader if better data becomes available

## What I Learned

**Critical Bug from My Previous Week:**
- I initially used placeholder A→G alleles for all variants → r≈0.03 (essentially zero)
- This scored wrong molecular changes (e.g., T→G instead of actual T→C variant)
- Spent a week debugging before discovering the issue!

**My Solution:**
- Fetch real alleles from dbSNP using rs IDs in QTL data
- Score actual variants → r=0.40 for caQTLs (13x improvement!)

**Other Key Insights:**
- Tissue matching is critical (CD4+ T cells for caQTLs, B cells for hQTLs)
- Use directional correlation, not absolute value (preserves effect sign)
- Deduplicate scorers when variants affect multiple related modalities
- Old datasets (pre-2015) often lack modern dbSNP annotations → can't fetch alleles

## Usage

```bash
# Full pipeline for both datasets
python scripts/01_prepare_qtls.py --datasets caQTLs hQTLs
python scripts/02_predict.py --datasets caQTLs hQTLs --force
python scripts/03_evaluate.py --datasets caQTLs hQTLs
python scripts/04_plots.py --datasets caQTLs hQTLs

# Test on small sample first
python scripts/01_prepare_qtls.py --datasets hQTLs --limit 100
python scripts/02_predict.py --datasets hQTLs --limit 100 --force
```
## Directory Structure
## Directory Structure

```
alphagenome-qtl-validation/
├── data/
│   ├── raw/               # Original QTL summary stats
│   ├── processed/         # Normalized QTLs with real alleles
│   └── genome/            # Reference genome (symlink)
├── scripts/
│   ├── 01_prepare_qtls.py    # Fetch real alleles from myvariant.info
│   ├── 02_predict.py          # Score variants with AlphaGenome
│   ├── 03_evaluate.py         # Compute correlations
│   └── 04_plots.py            # Generate visualizations
├── results/
│   ├── predictions/       # Variant predictions
│   ├── tables/            # Correlation metrics
│   └── plots/             # Scatter plots, distributions
└── logs/                  # Execution logs
```
