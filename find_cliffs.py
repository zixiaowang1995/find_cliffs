#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Ultra-Fast Activity Cliffs Detection Script with Command Line Interface
Based on MoleculeACE project optimization, specifically designed for classification tasks

Optimization strategies:
1. Early termination: Return immediately once similarity threshold is met
2. Pre-computed fingerprints: Avoid redundant molecular fingerprint calculations
3. Label-based grouping: Only compare compounds with different labels
4. Caching mechanism: Avoid redundant similarity calculations
5. Optimized data structures: Use more efficient data access patterns

Expected performance improvement: 2-5x speed boost
"""

import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, Scaffolds
from rdkit import DataStructs
from collections import defaultdict
import time
import argparse
from typing import Dict, List, Tuple, Set

class UltraFastCliffFinder:
    def __init__(self, similarity_threshold=0.9):
        self.similarity_threshold = similarity_threshold
        self.fingerprint_cache = {}  # Fingerprint cache
        self.scaffold_cache = {}     # Scaffold cache
        self.similarity_cache = {}   # Similarity cache
        
    def get_fingerprint(self, smiles: str):
        """Get molecular fingerprint with caching"""
        if smiles in self.fingerprint_cache:
            return self.fingerprint_cache[smiles]
            
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                self.fingerprint_cache[smiles] = None
                return None
            fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
            self.fingerprint_cache[smiles] = fp
            return fp
        except:
            self.fingerprint_cache[smiles] = None
            return None
    
    def get_scaffold(self, smiles: str):
        """Get molecular scaffold with caching"""
        if smiles in self.scaffold_cache:
            return self.scaffold_cache[smiles]
            
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                self.scaffold_cache[smiles] = None
                return None
            scaffold = Scaffolds.MurckoScaffold.GetScaffoldForMol(mol)
            if scaffold is None:
                self.scaffold_cache[smiles] = None
                return None
            scaffold_fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(scaffold, 2, nBits=2048)
            self.scaffold_cache[smiles] = scaffold_fp
            return scaffold_fp
        except:
            self.scaffold_cache[smiles] = None
            return None
    
    def calculate_tanimoto_similarity(self, fp1, fp2):
        """Calculate Tanimoto similarity"""
        if fp1 is None or fp2 is None:
            return 0.0
        return DataStructs.TanimotoSimilarity(fp1, fp2)
    
    def calculate_levenshtein_similarity(self, smiles1: str, smiles2: str):
        """Calculate Levenshtein similarity"""
        def levenshtein_distance(s1, s2):
            if len(s1) < len(s2):
                return levenshtein_distance(s2, s1)
            if len(s2) == 0:
                return len(s1)
            
            previous_row = list(range(len(s2) + 1))
            for i, c1 in enumerate(s1):
                current_row = [i + 1]
                for j, c2 in enumerate(s2):
                    insertions = previous_row[j + 1] + 1
                    deletions = current_row[j] + 1
                    substitutions = previous_row[j] + (c1 != c2)
                    current_row.append(min(insertions, deletions, substitutions))
                previous_row = current_row
            
            return previous_row[-1]
        
        max_len = max(len(smiles1), len(smiles2))
        if max_len == 0:
            return 1.0
        distance = levenshtein_distance(smiles1, smiles2)
        return 1.0 - (distance / max_len)
    
    def is_structurally_similar_fast(self, smiles1: str, smiles2: str):
        """Fast structural similarity check with early termination and caching"""
        # Check cache
        cache_key = tuple(sorted([smiles1, smiles2]))
        if cache_key in self.similarity_cache:
            return self.similarity_cache[cache_key]
        
        # 1. First check Tanimoto similarity (usually fastest)
        fp1 = self.get_fingerprint(smiles1)
        fp2 = self.get_fingerprint(smiles2)
        
        if fp1 is not None and fp2 is not None:
            tanimoto_sim = self.calculate_tanimoto_similarity(fp1, fp2)
            if tanimoto_sim >= self.similarity_threshold:
                self.similarity_cache[cache_key] = True
                return True
        
        # 2. Check scaffold similarity
        scaffold1 = self.get_scaffold(smiles1)
        scaffold2 = self.get_scaffold(smiles2)
        
        if scaffold1 is not None and scaffold2 is not None:
            scaffold_sim = self.calculate_tanimoto_similarity(scaffold1, scaffold2)
            if scaffold_sim >= self.similarity_threshold:
                self.similarity_cache[cache_key] = True
                return True
        
        # 3. Finally check Levenshtein similarity (slowest)
        levenshtein_sim = self.calculate_levenshtein_similarity(smiles1, smiles2)
        if levenshtein_sim >= self.similarity_threshold:
            self.similarity_cache[cache_key] = True
            return True
        
        self.similarity_cache[cache_key] = False
        return False
    
    def find_cliffs_optimized(self, df: pd.DataFrame, smiles_col: str, label_col: str):
        """Optimized activity cliffs finding algorithm"""
        print(f"Starting to process {len(df)} compounds...")
        
        # Group by labels, only compare compounds with different labels
        label_0_indices = df[df[label_col] == 0].index.tolist()
        label_1_indices = df[df[label_col] == 1].index.tolist()
        
        print(f"Label 0 compounds: {len(label_0_indices)}")
        print(f"Label 1 compounds: {len(label_1_indices)}")
        print(f"Total compound pairs to compare: {len(label_0_indices) * len(label_1_indices)}")
        
        cliff_pairs = set()  # Use set to avoid duplicates
        total_comparisons = len(label_0_indices) * len(label_1_indices)
        processed = 0
        
        # Pre-compute all required fingerprints
        print("Pre-computing molecular fingerprints...")
        all_smiles = set(df[smiles_col].tolist())
        for smiles in all_smiles:
            self.get_fingerprint(smiles)
            self.get_scaffold(smiles)
        
        print("Starting activity cliffs search...")
        start_time = time.time()
        
        # Only compare compounds with different labels
        for i, idx1 in enumerate(label_0_indices):
            smiles1 = df.loc[idx1, smiles_col]
            
            for j, idx2 in enumerate(label_1_indices):
                smiles2 = df.loc[idx2, smiles_col]
                
                # Skip identical SMILES
                if smiles1 == smiles2:
                    continue
                
                # Check structural similarity
                if self.is_structurally_similar_fast(smiles1, smiles2):
                    # Ensure smaller index comes first to avoid duplicates
                    pair = tuple(sorted([idx1, idx2]))
                    cliff_pairs.add(pair)
                
                processed += 1
                if processed % 10000 == 0:
                    elapsed = time.time() - start_time
                    progress = processed / total_comparisons * 100
                    print(f"Progress: {progress:.1f}% ({processed}/{total_comparisons}), "
                          f"Found {len(cliff_pairs)} activity cliff pairs, "
                          f"Time elapsed: {elapsed:.1f}s")
        
        elapsed = time.time() - start_time
        print(f"\nSearch completed!")
        print(f"Total time: {elapsed:.1f}s")
        print(f"Found {len(cliff_pairs)} activity cliff pairs")
        print(f"Cache statistics:")
        print(f"  Fingerprint cache: {len(self.fingerprint_cache)} entries")
        print(f"  Scaffold cache: {len(self.scaffold_cache)} entries")
        print(f"  Similarity cache: {len(self.similarity_cache)} entries")
        
        return cliff_pairs

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='Ultra-Fast Activity Cliffs Detection Tool',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument(
        '-i', '--input',
        type=str,
        default='1.csv',
        help='Input CSV file name'
    )
    
    parser.add_argument(
        '-o', '--output',
        type=str,
        default='activity_cliffs_ultra_fast.csv',
        help='Output CSV file name'
    )
    
    parser.add_argument(
        '-t', '--threshold',
        type=float,
        default=0.9,
        help='Similarity threshold for structural similarity'
    )
    
    parser.add_argument(
        '-s', '--smiles_column',
        type=str,
        default='smiles',
        help='Name of the SMILES column in the input file'
    )
    
    parser.add_argument(
        '-l', '--label_column',
        type=str,
        default='y',
        help='Name of the label column in the input file'
    )
    
    return parser.parse_args()

def main():
    # Parse command line arguments
    args = parse_arguments()
    
    print("=== Ultra-Fast Activity Cliffs Detection Tool ===")
    print(f"Input file: {args.input}")
    print(f"Output file: {args.output}")
    print(f"Similarity threshold: {args.threshold}")
    print(f"SMILES column: {args.smiles_column}")
    print(f"Label column: {args.label_column}")
    print()
    
    # Read data
    try:
        df = pd.read_csv(args.input)
        print(f"Successfully loaded data with {len(df)} rows")
    except Exception as e:
        print(f"Failed to read file: {e}")
        return
    
    # Check required columns
    if args.smiles_column not in df.columns:
        print(f"Error: SMILES column '{args.smiles_column}' not found")
        print(f"Available columns: {list(df.columns)}")
        return
    
    if args.label_column not in df.columns:
        print(f"Error: Label column '{args.label_column}' not found")
        print(f"Available columns: {list(df.columns)}")
        return
    
    # Check label values
    unique_labels = df[args.label_column].unique()
    print(f"Label values: {unique_labels}")
    
    if not (0 in unique_labels and 1 in unique_labels):
        print("Warning: Data does not contain both labels 0 and 1, cannot find activity cliffs")
        return
    
    # Create cliff finder and search for activity cliffs
    cliff_finder = UltraFastCliffFinder(similarity_threshold=args.threshold)
    cliff_pairs = cliff_finder.find_cliffs_optimized(df, args.smiles_column, args.label_column)
    
    if not cliff_pairs:
        print("No activity cliff compounds found")
        return
    
    # Prepare output data
    cliff_data = []
    for idx1, idx2 in cliff_pairs:
        row1 = df.loc[idx1].copy()
        row2 = df.loc[idx2].copy()
        
        # Add pairing information
        row1['pair_id'] = f"pair_{len(cliff_data)//2 + 1}"
        row1['pair_member'] = 'A'
        row1['partner_index'] = idx2
        
        row2['pair_id'] = f"pair_{len(cliff_data)//2 + 1}"
        row2['pair_member'] = 'B'
        row2['partner_index'] = idx1
        
        cliff_data.extend([row1, row2])
    
    # Create output DataFrame
    cliff_df = pd.DataFrame(cliff_data)
    
    # Save results
    try:
        cliff_df.to_csv(args.output, index=False)
        print(f"\nResults saved to {args.output}")
        print(f"Found {len(cliff_pairs)} activity cliff pairs, output {len(cliff_df)} rows of data")
        
        # Display statistics
        print("\n=== Result Statistics ===")
        print(f"Activity cliff pairs: {len(cliff_pairs)}")
        print(f"Total compounds involved: {len(cliff_df)}")
        print(f"Label 0 cliff compounds: {len(cliff_df[cliff_df[args.label_column] == 0])}")
        print(f"Label 1 cliff compounds: {len(cliff_df[cliff_df[args.label_column] == 1])}")
        
    except Exception as e:
        print(f"Failed to save file: {e}")

if __name__ == "__main__":
    main()