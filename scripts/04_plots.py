#!/usr/bin/env python3
"""
Step 4: Generate visualizations of AlphaGenome predictions vs QTL effects.

Creates:
1. Scatter plot with regression line (predicted vs observed)
2. Distribution plots (scores and betas)
3. Residual plots
4. Top variants analysis

Input:
- results/predictions/{dataset}.parquet

Output:
- results/plots/{dataset}_correlation.png
- results/plots/{dataset}_distributions.png
- results/plots/{dataset}_residuals.png
"""

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy.stats import spearmanr, pearsonr
from scipy import stats

# Set style
sns.set_style('whitegrid')
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.size'] = 10


def plot_correlation(df: pd.DataFrame, dataset: str, output_path: Path):
    """Create scatter plot of predictions vs QTL effects."""
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Scatter plot
    ax.scatter(df['quantile_score'], df['beta'], 
              alpha=0.5, s=20, color='steelblue', edgecolors='none')
    
    # Regression line
    slope, intercept, r_value, p_value, std_err = stats.linregress(df['quantile_score'], df['beta'])
    x_line = np.array([df['quantile_score'].min(), df['quantile_score'].max()])
    y_line = slope * x_line + intercept
    ax.plot(x_line, y_line, 'r-', linewidth=2, label=f'Linear fit (R²={r_value**2:.3f})')
    
    # Compute correlations
    spearman_r, spearman_p = spearmanr(df['quantile_score'], df['beta'])
    pearson_r, pearson_p = pearsonr(df['quantile_score'], df['beta'])
    
    # Labels and title
    ax.set_xlabel('AlphaGenome Quantile Score', fontsize=12)
    ax.set_ylabel('QTL Effect Size (Beta)', fontsize=12)
    ax.set_title(f'{dataset} - AlphaGenome Predictions vs QTL Effects\n'
                f'n={len(df):,} variants', fontsize=14, fontweight='bold')
    
    # Add stats box
    stats_text = (
        f'Spearman r = {spearman_r:.3f}\n'
        f'p = {spearman_p:.2e}\n\n'
        f'Pearson r = {pearson_r:.3f}\n'
        f'p = {pearson_p:.2e}'
    )
    ax.text(0.05, 0.95, stats_text, transform=ax.transAxes,
           verticalalignment='top', bbox=dict(boxstyle='round', 
           facecolor='wheat', alpha=0.8), fontsize=10)
    
    ax.legend(loc='lower right')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"  Saved correlation plot: {output_path}")


def plot_distributions(df: pd.DataFrame, dataset: str, output_path: Path):
    """Plot distributions of quantile scores and betas."""
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Quantile score distribution
    ax = axes[0, 0]
    ax.hist(df['quantile_score'], bins=50, color='steelblue', alpha=0.7, edgecolor='black')
    ax.axvline(df['quantile_score'].mean(), color='red', linestyle='--', 
              linewidth=2, label=f'Mean = {df["quantile_score"].mean():.3f}')
    ax.set_xlabel('Quantile Score', fontsize=11)
    ax.set_ylabel('Count', fontsize=11)
    ax.set_title('Distribution of AlphaGenome Scores', fontsize=12, fontweight='bold')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Beta distribution
    ax = axes[0, 1]
    ax.hist(df['beta'], bins=50, color='coral', alpha=0.7, edgecolor='black')
    ax.axvline(df['beta'].mean(), color='red', linestyle='--', 
              linewidth=2, label=f'Mean = {df["beta"].mean():.3f}')
    ax.set_xlabel('QTL Beta', fontsize=11)
    ax.set_ylabel('Count', fontsize=11)
    ax.set_title('Distribution of QTL Effect Sizes', fontsize=12, fontweight='bold')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Absolute values comparison
    ax = axes[1, 0]
    ax.hist(df['quantile_score'].abs(), bins=50, color='steelblue', 
           alpha=0.5, label='|Quantile Score|', edgecolor='black')
    ax.hist(df['beta'], bins=50, color='coral', alpha=0.5, 
           label='QTL Beta', edgecolor='black')
    ax.set_xlabel('Value', fontsize=11)
    ax.set_ylabel('Count', fontsize=11)
    ax.set_title('Magnitude Comparison', fontsize=12, fontweight='bold')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Sign agreement
    ax = axes[1, 1]
    df['same_sign'] = np.sign(df['quantile_score']) == np.sign(df['beta'])
    sign_counts = df['same_sign'].value_counts()
    colors = ['green' if x else 'red' for x in sign_counts.index]
    bars = ax.bar(['Same Sign', 'Opposite Sign'], 
                  [sign_counts.get(True, 0), sign_counts.get(False, 0)],
                  color=colors, alpha=0.7, edgecolor='black')
    ax.set_ylabel('Count', fontsize=11)
    ax.set_title(f'Directional Agreement\n{sign_counts.get(True, 0)/len(df)*100:.1f}% same sign', 
                fontsize=12, fontweight='bold')
    ax.grid(True, alpha=0.3, axis='y')
    
    # Add counts on bars
    for bar in bars:
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height,
               f'{int(height):,}', ha='center', va='bottom', fontsize=10)
    
    plt.suptitle(f'{dataset} - Distribution Analysis', fontsize=14, fontweight='bold', y=1.00)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"  Saved distribution plots: {output_path}")


def plot_residuals(df: pd.DataFrame, dataset: str, output_path: Path):
    """Plot residual analysis."""
    # Fit linear model
    slope, intercept, r_value, p_value, std_err = stats.linregress(df['quantile_score'], df['beta'])
    df['predicted_beta'] = slope * df['quantile_score'] + intercept
    df['residual'] = df['beta'] - df['predicted_beta']
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # Residual plot
    ax = axes[0]
    ax.scatter(df['predicted_beta'], df['residual'], 
              alpha=0.5, s=20, color='purple', edgecolors='none')
    ax.axhline(0, color='red', linestyle='--', linewidth=2)
    ax.set_xlabel('Predicted Beta', fontsize=11)
    ax.set_ylabel('Residual (Observed - Predicted)', fontsize=11)
    ax.set_title('Residual Plot', fontsize=12, fontweight='bold')
    ax.grid(True, alpha=0.3)
    
    # Residual distribution
    ax = axes[1]
    ax.hist(df['residual'], bins=50, color='purple', alpha=0.7, edgecolor='black')
    ax.axvline(0, color='red', linestyle='--', linewidth=2, label='Zero')
    ax.set_xlabel('Residual', fontsize=11)
    ax.set_ylabel('Count', fontsize=11)
    ax.set_title(f'Residual Distribution\nMean = {df["residual"].mean():.4f}, Std = {df["residual"].std():.4f}', 
                fontsize=12, fontweight='bold')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.suptitle(f'{dataset} - Residual Analysis', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"  Saved residual plots: {output_path}")


def plot_top_variants(df: pd.DataFrame, dataset: str, output_path: Path, n_top=20):
    """Plot top predicted and top observed variants."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 8))
    
    # Top by absolute quantile score
    ax = axes[0]
    df_sorted = df.nlargest(n_top, 'quantile_score', keep='all')[:n_top]
    colors = ['green' if x > 0 else 'red' for x in df_sorted['quantile_score']]
    y_pos = np.arange(len(df_sorted))
    ax.barh(y_pos, df_sorted['quantile_score'], color=colors, alpha=0.7, edgecolor='black')
    ax.set_yticks(y_pos)
    ax.set_yticklabels(df_sorted['variant_id'], fontsize=8)
    ax.set_xlabel('Quantile Score', fontsize=11)
    ax.set_title(f'Top {n_top} by AlphaGenome Score', fontsize=12, fontweight='bold')
    ax.grid(True, alpha=0.3, axis='x')
    ax.invert_yaxis()
    
    # Top by beta
    ax = axes[1]
    df_sorted = df.nlargest(n_top, 'beta', keep='all')[:n_top]
    y_pos = np.arange(len(df_sorted))
    ax.barh(y_pos, df_sorted['beta'], color='coral', alpha=0.7, edgecolor='black')
    ax.set_yticks(y_pos)
    ax.set_yticklabels(df_sorted['variant_id'], fontsize=8)
    ax.set_xlabel('QTL Beta', fontsize=11)
    ax.set_title(f'Top {n_top} by QTL Effect Size', fontsize=12, fontweight='bold')
    ax.grid(True, alpha=0.3, axis='x')
    ax.invert_yaxis()
    
    plt.suptitle(f'{dataset} - Top Variants', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"  Saved top variants plot: {output_path}")


def main():
    parser = argparse.ArgumentParser(description='Generate plots for AlphaGenome QTL benchmark')
    parser.add_argument('--datasets', nargs='+', default=['caQTLs'],
                       choices=['caQTLs', 'dsQTLs', 'eQTLs', 'hQTLs'])
    args = parser.parse_args()
    
    base_dir = Path(__file__).parent.parent
    predictions_dir = base_dir / "results/predictions"
    plots_dir = base_dir / "results/plots"
    plots_dir.mkdir(parents=True, exist_ok=True)
    
    for dataset in args.datasets:
        print(f"\n{'='*60}")
        print(f"Generating plots for {dataset}")
        print('='*60)
        
        # Load predictions
        pred_path = predictions_dir / f"{dataset}.parquet"
        if not pred_path.exists():
            print(f"  {pred_path} not found, skipping...")
            continue
        
        df = pd.read_parquet(pred_path)
        print(f"  Loaded {len(df):,} predictions")
        
        # Create plots
        plot_correlation(df, dataset, plots_dir / f"{dataset}_correlation.png")
        plot_distributions(df, dataset, plots_dir / f"{dataset}_distributions.png")
        plot_residuals(df, dataset, plots_dir / f"{dataset}_residuals.png")
        plot_top_variants(df, dataset, plots_dir / f"{dataset}_top_variants.png")
    
    print("\n✓ All plots generated")


if __name__ == '__main__':
    main()
