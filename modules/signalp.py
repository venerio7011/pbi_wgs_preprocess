#!/usr/bin/env python3
"""
Signal Peptide Prediction Module

Predicts signal peptides in bacterial protein sequences.
Simplified implementation for demonstration purposes.
"""

import os
from pathlib import Path
from typing import List, Dict, Tuple, Optional
import pandas as pd
import numpy as np
from tqdm import tqdm

from common import read_fasta, ensure_dir


class SignalPPredictor:
    """
    Signal Peptide Predictor

    Predicts presence and cleavage site of signal peptides
    using hydrophobicity and charge patterns.
    """

    # Kyte-Doolittle hydrophobicity scale
    HYDROPHOBICITY = {
        'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
        'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
        'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
        'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2,
    }

    # Charged amino acids
    POSITIVE_AA = set('RKH')
    NEGATIVE_AA = set('DE')

    def __init__(self, min_length: int = 15, max_length: int = 30):
        """
        Initialize Signal Peptide Predictor.

        Args:
            min_length: Minimum signal peptide length
            max_length: Maximum signal peptide length
        """
        self.min_length = min_length
        self.max_length = max_length

    def _calculate_hydrophobicity(self, sequence: str) -> float:
        """Calculate average hydrophobicity of a sequence."""
        if not sequence:
            return 0.0
        seq_upper = sequence.upper()
        return sum(self.HYDROPHOBICITY.get(aa, 0) for aa in seq_upper) / len(seq_upper)

    def _count_charged(self, sequence: str) -> Tuple[int, int]:
        """Count positive and negative charged residues."""
        seq_upper = sequence.upper()
        positive = sum(1 for aa in seq_upper if aa in self.POSITIVE_AA)
        negative = sum(1 for aa in seq_upper if aa in self.NEGATIVE_AA)
        return positive, negative

    def _find_cleavage_site(self, sequence: str) -> Tuple[int, float]:
        """
        Find potential signal peptide cleavage site.

        Returns:
            Tuple of (cleavage_position, score)
        """
        best_pos = 0
        best_score = 0.0

        seq_upper = sequence.upper()
        search_end = min(len(seq_upper), self.max_length + 5)

        for pos in range(self.min_length, search_end):
            # Check for Ala-X-Ala motif (common cleavage pattern)
            if pos >= 3:
                if seq_upper[pos-3] == 'A' and seq_upper[pos-1] == 'A':
                    score = 0.8
                elif seq_upper[pos-1] in 'AGS':
                    score = 0.5
                else:
                    score = 0.2

                # N-region should be positively charged
                n_region = seq_upper[:8]
                pos_count, _ = self._count_charged(n_region)
                if pos_count >= 1:
                    score += 0.2

                # H-region should be hydrophobic
                h_start = min(8, pos - 10)
                h_end = min(pos - 3, pos)
                if h_end > h_start:
                    h_region = seq_upper[h_start:h_end]
                    h_hydro = self._calculate_hydrophobicity(h_region)
                    if h_hydro > 1.0:
                        score += 0.3

                if score > best_score:
                    best_score = score
                    best_pos = pos

        return best_pos, best_score

    def predict(
        self,
        sequences: List[Tuple[str, str]],
        threshold: float = 0.5,
    ) -> pd.DataFrame:
        """
        Predict signal peptides for sequences.

        Args:
            sequences: List of (header, sequence) tuples
            threshold: Score threshold for positive prediction

        Returns:
            DataFrame with predictions
        """
        print(f"[SignalP] Predicting signal peptides for {len(sequences)} sequences...")

        results = []
        for header, seq in tqdm(sequences, desc="Analyzing sequences"):
            protein_id = header.split()[0]

            # Find cleavage site
            cleavage_pos, score = self._find_cleavage_site(seq)

            # Determine prediction
            has_sp = score >= threshold

            results.append({
                'protein_id': protein_id,
                'description': header,
                'length': len(seq),
                'has_signal_peptide': has_sp,
                'cleavage_site': cleavage_pos if has_sp else 0,
                'sp_score': score,
                'sp_sequence': seq[:cleavage_pos] if has_sp else '',
            })

        return pd.DataFrame(results)

    def run(
        self,
        input_file: str,
        output_dir: str,
    ) -> pd.DataFrame:
        """
        Run signal peptide prediction pipeline.

        Args:
            input_file: Input FASTA file
            output_dir: Output directory

        Returns:
            DataFrame with predictions
        """
        ensure_dir(output_dir)

        sequences = read_fasta(input_file)
        print(f"[SignalP] Loaded {len(sequences)} sequences from {input_file}")

        results = self.predict(sequences)

        output_file = Path(output_dir) / "signalp_predictions.tsv"
        results.to_csv(output_file, sep='\t', index=False)
        print(f"[SignalP] Results saved to {output_file}")

        sp_count = results['has_signal_peptide'].sum()
        print(f"\n[SignalP] Proteins with signal peptide: {sp_count}/{len(results)}")

        return results
