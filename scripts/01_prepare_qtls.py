#!/usr/bin/env python3
"""
Step 1: Prepare QTL data with REAL alleles from dbSNP/myvariant.info

Fetches actual ref/alt alleles for each QTL variant using rs IDs.
This fixes the root cause: previous approach used placeholder A→G for all variants.

Input:
- data/raw/caQTLs_GSE86886/ATAC-QTLs.csv (has rs IDs)
- data/raw/hQTLs_GSE116193/Pelikan_et_al_hQTL_summary.csv

Output:
- data/processed/{dataset}.parquet
  Columns: variant_id, chrom, pos, ref, alt, beta, modality
"""

import pandas as pd
import requests
import time
from pathlib import Path
from tqdm import tqdm
import argparse


def extract_rs_id(snp_str: str) -> str:
    """Extract rs ID from SNP string like 'chr12:9436083_rs61916194'."""
    if '_rs' in str(snp_str):
        return 'rs' + str(snp_str).split('_rs')[1]
    return None


def fetch_alleles_batch(rs_ids: list, assembly='hg38') -> dict:
    """
    Fetch alleles for multiple rs IDs from myvariant.info API.
    
    Returns: {rsid: {'ref': 'C', 'alt': 'T', 'chrom': '12', 'pos': 9283487}}
    """
    url = 'https://myvariant.info/v1/variant'
    ids_str = ','.join(rs_ids)
    
    try:
        response = requests.post(
            url,
            data={'ids': ids_str, 'assembly': assembly, 'fields': 'dbsnp.ref,dbsnp.alt,chrom,dbsnp.hg38.start'},
            timeout=30
        )
        response.raise_for_status()
        data = response.json()
        
        results = {}
        for item in data:
            if 'dbsnp' in item and 'query' in item:
                rsid = item['query']
                dbsnp = item['dbsnp']
                ref = dbsnp.get('ref')
                alt = dbsnp.get('alt')
                
                if ref and alt and isinstance(ref, str) and isinstance(alt, str):
                    results[rsid] = {
                        'ref': ref,
                        'alt': alt,
                        'chrom': item.get('chrom'),
                        'pos': dbsnp.get('hg38', {}).get('start')
                    }
        return results
    except Exception as e:
        print(f"  Error fetching batch: {e}")
        return {}


def load_caQTLs(base_dir: Path) -> pd.DataFrame:
    """Load caQTLs with real alleles from myvariant.info."""
    csv_path = base_dir / "data/raw/caQTLs_GSE86886/ATAC-QTLs.csv"
    df = pd.read_csv(csv_path)
    
    print(f"Loading {len(df)} caQTLs...")
    
    # Extract rs IDs
    rs_ids = []
    rs_id_map = {}
    
    for idx, row in df.iterrows():
        snp = str(row['SNP'])
        rs_id = extract_rs_id(snp)
        if rs_id:
            variant_id = f"caQTL_{idx}"
            rs_ids.append(rs_id)
            rs_id_map[variant_id] = rs_id
    
    print(f"Found {len(rs_ids)} rs IDs, fetching alleles...")
    
    # Fetch alleles in batches
    batch_size = 100
    all_alleles = {}
    
    for i in tqdm(range(0, len(rs_ids), batch_size), desc="Fetching alleles"):
        batch = rs_ids[i:i+batch_size]
        alleles = fetch_alleles_batch(batch)
        all_alleles.update(alleles)
        time.sleep(0.5)  # Rate limiting
    
    print(f"Successfully fetched alleles for {len(all_alleles)}/{len(rs_ids)} variants")
    
    # Build records with real alleles
    records = []
    for idx, row in df.iterrows():
        variant_id = f"caQTL_{idx}"
        rs_id = rs_id_map.get(variant_id)
        
        # Try to get real alleles from API
        if rs_id and rs_id in all_alleles:
            allele_info = all_alleles[rs_id]
            chrom = allele_info['chrom']
            pos = allele_info['pos']
            ref = allele_info['ref']
            alt = allele_info['alt']
        else:
            # Fallback: parse from SNP field
            snp = str(row['SNP'])
            if ':' in snp:
                parts = snp.split(':')
                chrom = parts[0].replace('chr', '')
                try:
                    pos = int(parts[1].split('_')[0])
                except:
                    continue
                ref, alt = 'N', 'N'  # Unknown
            else:
                continue
        
        records.append({
            'variant_id': variant_id,
            'chrom': chrom,
            'pos': pos,
            'ref': ref,
            'alt': alt,
            'beta': float(row.get('beta', row.get('Beta', 0))),
            'modality': 'ATAC,DNase,H3K27ac'
        })
    
    return pd.DataFrame(records)


def load_hQTLs(base_dir: Path, limit: int = None) -> pd.DataFrame:
    """Load hQTLs with real alleles (already in the file + fetch from API)."""
    csv_path = base_dir / "data/raw/hQTLs_GSE116193/Pelikan_et_al_hQTL_summary.csv"
    df = pd.read_csv(csv_path)
    
    if limit:
        df = df.head(limit)
        print(f"Testing on first {len(df)} hQTLs...")
    else:
        print(f"Loading {len(df)} hQTLs...")
    
    # hQTLs have ref/alt alleles already, but positions are hg19
    # We'll use rs IDs to get hg38 positions from myvariant.info
    
    # Extract rs IDs
    rs_ids = []
    rs_id_map = {}
    
    for idx, row in df.iterrows():
        rs_id = str(row['epiQTL rsID']).strip()
        if rs_id.startswith('rs'):
            variant_id = f"hQTL_{idx}"
            rs_ids.append(rs_id)
            rs_id_map[variant_id] = {
                'rs_id': rs_id,
                'chrom': str(row['Chr']),
                'pos_hg19': int(row['Bp (hg19)']),
                'ref': str(row['Ref Allelec']).strip(),
                'alt': str(row['Alt Alleled']).strip(),
                'effect_size': float(row['log2 (effect size)']),
                'has_H3K27ac': pd.notna(row['H3K27ac p-valueb']) and row['H3K27ac p-valueb'] != '.',
                'has_H3K4me1': pd.notna(row['H3K4me1 p-valueb']) and row['H3K4me1 p-valueb'] != '.'
            }
    
    print(f"Found {len(rs_ids)} variants with rs IDs, fetching hg38 positions...")
    
    # Fetch hg38 positions in batches
    batch_size = 100
    all_positions = {}
    
    for i in tqdm(range(0, len(rs_ids), batch_size), desc="Fetching hg38 positions"):
        batch = rs_ids[i:i+batch_size]
        alleles = fetch_alleles_batch(batch, assembly='hg38')
        
        for rs_id, info in alleles.items():
            if 'pos' in info and info['pos']:
                all_positions[rs_id] = {
                    'chrom': info.get('chrom'),
                    'pos_hg38': info['pos']
                }
        
        time.sleep(0.5)
    
    print(f"Successfully fetched hg38 positions for {len(all_positions)}/{len(rs_ids)} variants ({len(all_positions)/len(rs_ids)*100:.1f}%)")
    
    # Build records
    records = []
    for variant_id, info in rs_id_map.items():
        rs_id = info['rs_id']
        
        # Use hg38 position if available, otherwise keep hg19
        if rs_id in all_positions:
            chrom = all_positions[rs_id]['chrom']
            pos = all_positions[rs_id]['pos_hg38']
        else:
            chrom = info['chrom']
            pos = info['pos_hg19']
        
        # Determine which modality to use
        modalities = []
        if info['has_H3K27ac']:
            modalities.append('H3K27ac')
        if info['has_H3K4me1']:
            modalities.append('H3K4me1')
        
        if not modalities:
            continue  # Skip if no significant QTL
        
        records.append({
            'variant_id': variant_id,
            'chrom': chrom,
            'pos': pos,
            'ref': info['ref'],
            'alt': info['alt'],
            'beta': info['effect_size'],  # Directional log2 fold change (preserve sign)
            'modality': ','.join(modalities)
        })
    
    return pd.DataFrame(records)


def load_eQTLs(base_dir: Path) -> pd.DataFrame:
    """Load eQTLs (similar pattern to caQTLs)."""
    csv_path = base_dir / "data/raw/caQTLs_GSE86886/eQTLs.csv"
    
    if not csv_path.exists():
        return pd.DataFrame()
    
    df = pd.read_csv(csv_path)
    print(f"Loading {len(df)} eQTLs...")
    
    # Similar rs ID extraction and fetching logic as caQTLs
    # For brevity, using simpler approach here
    records = []
    for idx, row in df.iterrows():
        snp = str(row.get('SNP', ''))
        if ':' in snp:
            parts = snp.split(':')
            chrom = parts[0].replace('chr', '')
            try:
                pos = int(parts[1].split('_')[0])
            except:
                continue
                
            records.append({
                'variant_id': f"eQTL_{idx}",
                'chrom': chrom,
                'pos': pos,
                'ref': 'N',  # Would need similar API lookup
                'alt': 'N',
                'beta': float(row.get('beta', 0)),
                'modality': 'RNA_SEQ'
            })
    
    return pd.DataFrame(records)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--datasets', nargs='+', default=['caQTLs'],
                       choices=['caQTLs', 'hQTLs'])
    parser.add_argument('--limit', type=int, default=None,
                       help='Test on first N variants only')
    args = parser.parse_args()
    
    base_dir = Path(__file__).parent.parent
    output_dir = base_dir / "data/processed"
    output_dir.mkdir(exist_ok=True)
    
    for dataset in args.datasets:
        print(f"\n{'='*60}")
        print(f"Processing {dataset}")
        print('='*60)
        
        if dataset == 'caQTLs':
            df = load_caQTLs(base_dir)
        elif dataset == 'hQTLs':
            df = load_hQTLs(base_dir, limit=args.limit)
        else:
            print(f"  {dataset} not yet implemented")
            continue
        
        if len(df) == 0:
            print(f"  No data for {dataset}")
            continue
        
        # Save
        output_path = output_dir / f"{dataset}.parquet"
        df.to_parquet(output_path, index=False)
        
        print(f"\n  Saved {len(df)} variants to {output_path}")
        print(f"  Real alleles: {(df['ref'] != 'N').sum()}/{len(df)}")
        print(f"\n  Sample:")
        print(df[['variant_id', 'chrom', 'pos', 'ref', 'alt', 'beta']].head())
    
    print("\n✓ Data preparation complete")


if __name__ == '__main__':
    main()
