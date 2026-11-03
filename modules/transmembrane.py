#!/usr/bin/env python3
"""
Transmembrane Topology Prediction Module

Predicts transmembrane regions in bacterial proteins.
Simplified implementation for demonstration purposes.
"""

import os
from pathlib import Path
from typing import List, Dict, Tuple, Optional
import pandas as pd
import numpy as np
from tqdm import tqdm

from common import read_fasta, ensure_dir


class TransmembranePredictor:
    """
    Transmembrane Topology Predictor

    Predicts:
    - Number of transmembrane helices
    - Topology (in/out orientation)
    - Transmembrane regions
    """

    # Kyte-Doolittle hydrophobicity scale
    HYDROPHOBICITY = {
        'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
        'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
        'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
        'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2,
    }

    def __init__(self, window_size: int = 19, threshold: float = 1.6):
        """
        Initialize Transmembrane Predictor.

        Args:
            window_size: Window size for hydrophobicity calculation
            threshold: Hydrophobicity threshold for TM region
        """
        self.window_size = window_size
        self.threshold = threshold

    def _calculate_hydrophobicity_profile(self, sequence: str) -> List[float]:
        """Calculate hydrophobicity profile using sliding window."""
        seq_upper = sequence.upper()
        profile = []
        half_window = self.window_size // 2

        for i in range(len(seq_upper)):
            start = max(0, i - half_window)
            end = min(len(seq_upper), i + half_window + 1)
            window = seq_upper[start:end]

            if len(window) > 0:
                hydro = sum(self.HYDROPHOBICITY.get(aa, 0) for aa in window) / len(window)
            else:
                hydro = 0

            profile.append(hydro)

        return profile

    def _find_tm_regions(self, profile: List[float], min_length: int = 15) -> List[Tuple[int, int]]:
        """Find transmembrane regions from hydrophobicity profile."""
        regions = []
        in_region = False
        start = 0

        for i, hydro in enumerate(profile):
            if hydro >= self.threshold:
                if not in_region:
                    in_region = True
                    start = i
            else:
                if in_region:
                    if i - start >= min_length:
                        regions.append((start, i))
                    in_region = False

        if in_region and len(profile) - start >= min_length:
            regions.append((start, len(profile)))

        return regions

    def _classify_protein(self, tm_count: int, sequence: str) -> str:
        """Classify protein based on TM helix count."""
        if tm_count == 0:
            seq_upper = sequence.upper()
            aromatic_ratio = sum(1 for aa in seq_upper if aa in 'FWY') / len(seq_upper)
            if aromatic_ratio > 0.1 and len(seq_upper) > 200:
                return 'BETA'
            return 'GLOBULAR'
        elif tm_count == 1:
            return 'TM1'
        elif tm_count <= 6:
            return f'TM{tm_count}'
        else:
            return 'TM7+'

    def predict(
        self,
        sequences: List[Tuple[str, str]],
    ) -> pd.DataFrame:
        """
        Predict transmembrane topology for sequences.

        Args:
            sequences: List of (header, sequence) tuples

        Returns:
            DataFrame with predictions
        """
        print(f"[TM] Predicting transmembrane topology for {len(sequences)} sequences...")

        results = []
        for header, seq in tqdm(sequences, desc="Analyzing topology"):
            protein_id = header.split()[0]

            profile = self._calculate_hydrophobicity_profile(seq)
            tm_regions = self._find_tm_regions(profile)
            tm_count = len(tm_regions)
            classification = self._classify_protein(tm_count, seq)
            tm_regions_str = ';'.join(f"{s}-{e}" for s, e in tm_regions) if tm_regions else ''
            avg_hydro = np.mean(profile) if profile else 0

            results.append({
                'protein_id': protein_id,
                'description': header,
                'length': len(seq),
                'tm_count': tm_count,
                'classification': classification,
                'tm_regions': tm_regions_str,
                'avg_hydrophobicity': avg_hydro,
                'is_membrane_protein': tm_count > 0 or classification == 'BETA',
            })

        return pd.DataFrame(results)

    def run(
        self,
        input_file: str,
        output_dir: str,
    ) -> pd.DataFrame:
        """
        Run transmembrane prediction pipeline.

        Args:
            input_file: Input FASTA file
            output_dir: Output directory

        Returns:
            DataFrame with predictions
        """
        ensure_dir(output_dir)

        sequences = read_fasta(input_file)
        print(f"[TM] Loaded {len(sequences)} sequences from {input_file}")

        results = self.predict(sequences)

        output_file = Path(output_dir) / "transmembrane_predictions.tsv"
        results.to_csv(output_file, sep='\t', index=False)
        print(f"[TM] Results saved to {output_file}")

        print(f"\n[TM] Classification summary:")
        for cls in results['classification'].unique():
            count = (results['classification'] == cls).sum()
            print(f"  {cls}: {count}")

        membrane_count = results['is_membrane_protein'].sum()
        print(f"\n[TM] Membrane proteins: {membrane_count}/{len(results)}")

        return results
