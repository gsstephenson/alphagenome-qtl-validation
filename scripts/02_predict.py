#!/usr/bin/env python3
"""
Step 2: Run AlphaGenome variant scoring with REAL alleles and CD4+ T cell filtering.

Following the exact pattern from alphagenome_examples/batch_variant_scoring.ipynb

Input:
- data/processed/{dataset}.parquet (with real ref/alt alleles)

Output:
- results/predictions/{dataset}.parquet
  Columns: variant_id, chrom, pos, ref, alt, modality, quantile_score, beta
"""

import argparse
import pandas as pd
import numpy as np
import os
from pathlib import Path
from tqdm import tqdm


def load_alphagenome_client(api_key: str):
    """Initialize AlphaGenome client."""
    from alphagenome.data import genome
    from alphagenome.models import dna_client, variant_scorers
    
    client = dna_client.create(api_key=api_key)
    return client, genome, variant_scorers, dna_client


def get_modality_scorers(variant_scorers, modalities: str):
    """Get RECOMMENDED_VARIANT_SCORERS for specified modalities."""
    modality_map = {
        'ATAC': 'ATAC',
        'DNase': 'DNASE',
        'H3K27ac': 'CHIP_HISTONE',
        'H3K4me1': 'CHIP_HISTONE',
        'RNA_SEQ': 'RNA_SEQ'
    }
    
    # Use set to deduplicate scorer keys (e.g., H3K27ac and H3K4me1 both -> CHIP_HISTONE)
    scorer_keys = set()
    for mod in modalities.split(','):
        mod = mod.strip()
        if mod in modality_map:
            scorer_keys.add(modality_map[mod])
    
    # Get unique scorers
    scorers = [variant_scorers.RECOMMENDED_VARIANT_SCORERS[key] for key in scorer_keys]
    
    return scorers


def score_variants_batch(client, genome, variant_scorers_module, dna_client,
                         variants_df: pd.DataFrame, modalities: str, tissue_curie: str = None) -> pd.DataFrame:
    """
    Score variants using AlphaGenome following official batch_variant_scoring pattern.
    
    Key differences from old approach:
    1. Using REAL ref/alt alleles (not placeholders)
    2. Proper tissue filtering (CD4+ T cells for caQTLs, B cells for hQTLs)
    3. Aggregating across tracks per variant-modality pair
    """
    # Tissue ontology curies
    TISSUES = {
        'CD4_T_CELL': 'CL:0000624',
        'B_CELL': 'CL:0000236'
    }
    
    # Default to CD4+ T cells if not specified
    if tissue_curie is None:
        tissue_curie = TISSUES['CD4_T_CELL']
    
    # Get scorers for requested modalities
    scorers = get_modality_scorers(variant_scorers_module, modalities)
    
    if not scorers:
        raise ValueError(f"No valid scorers for modalities: {modalities}")
    
    print(f"  Using {len(scorers)} scorers for: {modalities}")
    
    # Score each variant (following batch_variant_scoring.ipynb pattern)
    results = []
    
    for _, row in tqdm(variants_df.iterrows(), total=len(variants_df), desc="  Scoring variants"):
        # Skip if alleles unknown
        if row['ref'] == 'N' or row['alt'] == 'N':
            continue
        
        try:
            # Create variant
            variant = genome.Variant(
                chromosome=f"chr{row['chrom']}" if not str(row['chrom']).startswith('chr') else str(row['chrom']),
                position=int(row['pos']),
                reference_bases=row['ref'],
                alternate_bases=row['alt'],
                name=row['variant_id']
            )
            
            # Create 1MB interval
            interval = variant.reference_interval.resize(dna_client.SEQUENCE_LENGTH_1MB)
            
            # Score variant
            variant_scores = client.score_variant(
                interval=interval,
                variant=variant,
                variant_scorers=scorers,
                organism=dna_client.Organism.HOMO_SAPIENS
            )
            
            results.append(variant_scores)
            
        except Exception as e:
            print(f"    Warning: Failed {row['variant_id']}: {e}")
            continue
    
    if not results:
        raise ValueError("No variants successfully scored!")
    
    # Tidy scores (converts to per-tissue DataFrame)
    print(f"  Converting {len(results)} results to tidy format...")
    df = variant_scorers_module.tidy_scores(results)
    
    # Filter to specified tissue
    tissue_name = [k for k, v in TISSUES.items() if v == tissue_curie][0]
    print(f"  Filtering to {tissue_name} ({tissue_curie}) (before: {len(df)} rows)...")
    df = df[df['ontology_curie'] == tissue_curie].copy()
    print(f"  After filtering: {len(df)} rows from {df['biosample_name'].nunique()} tracks")
    
    if len(df) == 0:
        raise ValueError(f"No {tissue_name} predictions found!")
    
    # Extract variant_id as string (it's currently a Variant object)
    df['variant_id_str'] = df['variant_id'].astype(str)
    
    # Aggregate across tracks per variant
    # Group by variant_id and take mean quantile_score across CD4+ tracks
    agg_df = df.groupby(['variant_id_str'], as_index=False).agg({
        'quantile_score': 'mean',
        'raw_score': 'mean'
    })
    
    # Merge back with variant info (beta, modality)
    # The variant_id_str format is "chr12:9283487:C>T"
    # We need to match it to our variant_id like "caQTL_0"
    # Create a lookup dict from position to variant_id
    variant_lookup = {}
    for _, row in variants_df.iterrows():
        # Skip if position is NaN
        if pd.isna(row['pos']):
            continue
        key = f"chr{row['chrom']}:{int(row['pos'])}:{row['ref']}>{row['alt']}"
        variant_lookup[key] = {
            'variant_id': row['variant_id'],
            'chrom': row['chrom'],
            'pos': row['pos'],
            'ref': row['ref'],
            'alt': row['alt'],
            'beta': row['beta'],
            'modality': row['modality']
        }
    
    # Map back to original variant IDs
    for idx, row in agg_df.iterrows():
        key = row['variant_id_str']
        if key in variant_lookup:
            for col, val in variant_lookup[key].items():
                agg_df.at[idx, col] = val
    
    # Drop the temp column and remove unmapped rows
    agg_df = agg_df.drop(columns=['variant_id_str'])
    agg_df = agg_df.dropna(subset=['variant_id'])  # Remove rows without mapping
    
    print(f"  Successfully mapped {len(agg_df)} variants")
    
    return agg_df


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--datasets', nargs='+', default=['caQTLs'],
                       choices=['caQTLs', 'hQTLs'])
    parser.add_argument('--limit', type=int, help='Test on first N variants')
    parser.add_argument('--force', action='store_true')
    args = parser.parse_args()
    
    # Get API key
    api_key = os.getenv('ALPHA_GENOME_KEY')
    if not api_key:
        raise ValueError("ALPHA_GENOME_KEY environment variable not set")
    
    base_dir = Path(__file__).parent.parent
    output_dir = base_dir / "results/predictions"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Initialize AlphaGenome
    print("Initializing AlphaGenome...")
    client, genome, variant_scorers, dna_client = load_alphagenome_client(api_key)
    
    for dataset in args.datasets:
        print(f"\n{'='*60}")
        print(f"Scoring {dataset}")
        print('='*60)
        
        # Load processed variants (with real alleles!)
        input_path = base_dir / f"data/processed/{dataset}.parquet"
        if not input_path.exists():
            print(f"  {input_path} not found, run 01_prepare_qtls.py first")
            continue
        
        output_path = output_dir / f"{dataset}.parquet"
        if output_path.exists() and not args.force:
            print(f"  {output_path} exists, use --force to overwrite")
            continue
        
        variants_df = pd.read_parquet(input_path)
        
        if args.limit:
            variants_df = variants_df.head(args.limit)
        
        print(f"  Loaded {len(variants_df)} variants")
        print(f"  Real alleles: {(variants_df['ref'] != 'N').sum()}/{len(variants_df)}")
        
        # Determine tissue curie based on dataset
        # caQTLs: CD4+ T cells (CL:0000624) - Nedelec et al.
        # hQTLs: B cells (CL:0000236) - Pelikan et al. LCLs
        tissue_curie = 'CL:0000624' if dataset == 'caQTLs' else 'CL:0000236'
        
        # Group variants by modality (e.g., hQTLs have H3K27ac, H3K4me1, or both)
        modality_groups = variants_df.groupby('modality')
        print(f"  Found {len(modality_groups)} modality groups:")
        for modality, group_df in modality_groups:
            print(f"    - {modality}: {len(group_df)} variants")
        
        # Score each modality group separately
        all_predictions = []
        for modality, group_df in modality_groups:
            print(f"\n  Scoring {len(group_df)} variants with modality: {modality}")
            group_predictions = score_variants_batch(
                client, genome, variant_scorers, dna_client,
                group_df, modality, tissue_curie
            )
            all_predictions.append(group_predictions)
        
        # Combine predictions from all modality groups
        predictions_df = pd.concat(all_predictions, ignore_index=True)
        
        # Save
        predictions_df.to_parquet(output_path, index=False)
        
        print(f"\n  Saved {len(predictions_df)} predictions to {output_path}")
        print(f"  Score range: [{predictions_df['quantile_score'].min():.3f}, "
              f"{predictions_df['quantile_score'].max():.3f}]")
    
    print("\n✓ Prediction complete")


if __name__ == '__main__':
    main()
