#!/usr/bin/env python3
"""
Step 3: Evaluate AlphaGenome predictions against QTL effect sizes.

Computes correlation between AlphaGenome quantile_score and QTL beta.

Input:
- results/predictions/{dataset}.parquet

Output:
- results/tables/{dataset}_metrics.txt
"""

import argparse
import pandas as pd
from pathlib import Path
from scipy.stats import spearmanr, pearsonr
import numpy as np


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--datasets', nargs='+', default=['caQTLs'])
    args = parser.parse_args()
    
    base_dir = Path(__file__).parent.parent
    output_dir = base_dir / "results/tables"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    for dataset in args.datasets:
        pred_path = base_dir / f"results/predictions/{dataset}.parquet"
        
        if not pred_path.exists():
            print(f"{dataset}: predictions not found")
            continue
        
        df = pd.read_parquet(pred_path)
        
        # Use raw quantile_score (directional) vs beta
        # Both are directional effect sizes
        
        # Compute correlations
        spearman_r, spearman_p = spearmanr(df['quantile_score'], df['beta'])
        pearson_r, pearson_p = pearsonr(df['quantile_score'], df['beta'])
        
        # Print results
        print(f"\n{dataset} Results (n={len(df)}):")
        print(f"  Spearman r: {spearman_r:.4f} (p={spearman_p:.2e})")
        print(f"  Pearson r:  {pearson_r:.4f} (p={pearson_p:.2e})")
        print(f"  Note: Correlating quantile_score with beta (directional effects)")
        
        # Save
        output_path = output_dir / f"{dataset}_metrics.txt"
        with open(output_path, 'w') as f:
            f.write(f"{dataset} Evaluation\n")
            f.write(f"{'='*60}\n")
            f.write(f"n_variants: {len(df)}\n")
            f.write(f"spearman_r: {spearman_r:.4f}\n")
            f.write(f"spearman_p: {spearman_p:.2e}\n")
            f.write(f"pearson_r: {pearson_r:.4f}\n")
            f.write(f"pearson_p: {pearson_p:.2e}\n")
            f.write(f"\nNote: Correlating quantile_score with beta (directional effects)\n")
            f.write(f"\nScore stats:\n")
            f.write(f"  mean: {df['quantile_score'].mean():.4f}\n")
            f.write(f"  std: {df['quantile_score'].std():.4f}\n")
            f.write(f"  range: [{df['quantile_score'].min():.4f}, {df['quantile_score'].max():.4f}]\n")
            f.write(f"  abs_mean: {df['quantile_score'].abs().mean():.4f}\n")
        
        print(f"  Saved to {output_path}")


if __name__ == '__main__':
    main()
