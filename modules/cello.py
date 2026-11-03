#!/usr/bin/env python3
"""
CELLO-like Subcellular Localization Prediction Module

Predicts subcellular localization of bacterial proteins.
Simplified implementation for demonstration purposes.
"""

import os
from pathlib import Path
from typing import List, Dict, Tuple, Optional
import pandas as pd
import numpy as np
from tqdm import tqdm

from common import read_fasta, ensure_dir


class CELLOPredictor:
    """
    CELLO-like Subcellular Localization Predictor

    Predicts protein localization based on amino acid composition
    and sequence features.
    """

    # Localization categories for Gram-negative bacteria
    LOCALIZATIONS = [
        'Cytoplasmic',
        'InnerMembrane',
        'Periplasmic',
        'OuterMembrane',
        'Extracellular',
    ]

    # Feature weights for each localization
    FEATURE_WEIGHTS = {
        'Cytoplasmic': {'hydrophobic': -0.3, 'charged': 0.4, 'aromatic': 0.1},
        'InnerMembrane': {'hydrophobic': 0.6, 'charged': -0.2, 'aromatic': 0.2},
        'Periplasmic': {'hydrophobic': -0.2, 'charged': 0.3, 'aromatic': 0.1},
        'OuterMembrane': {'hydrophobic': 0.4, 'charged': -0.1, 'aromatic': 0.5},
        'Extracellular': {'hydrophobic': -0.4, 'charged': 0.2, 'aromatic': 0.0},
    }

    def __init__(self):
        """Initialize CELLO Predictor."""
        pass

    def _calculate_aa_composition(self, sequence: str) -> Dict[str, float]:
        """Calculate amino acid composition features."""
        seq_upper = sequence.upper()
        length = len(seq_upper)

        if length == 0:
            return {'hydrophobic': 0, 'charged': 0, 'aromatic': 0}

        hydrophobic = sum(1 for aa in seq_upper if aa in 'AILMFVPW') / length
        charged = sum(1 for aa in seq_upper if aa in 'RKDEH') / length
        aromatic = sum(1 for aa in seq_upper if aa in 'FWY') / length

        return {
            'hydrophobic': hydrophobic,
            'charged': charged,
            'aromatic': aromatic,
        }

    def _calculate_scores(self, features: Dict[str, float]) -> Dict[str, float]:
        """Calculate localization scores based on features."""
        scores = {}

        for loc, weights in self.FEATURE_WEIGHTS.items():
            score = 0.5
            for feat, weight in weights.items():
                score += features.get(feat, 0) * weight
            scores[loc] = max(0, min(1, score))

        # Normalize scores
        total = sum(scores.values())
        if total > 0:
            scores = {k: v / total for k, v in scores.items()}

        return scores

    def predict(
        self,
        sequences: List[Tuple[str, str]],
    ) -> pd.DataFrame:
        """
        Predict subcellular localization for sequences.

        Args:
            sequences: List of (header, sequence) tuples

        Returns:
            DataFrame with predictions
        """
        print(f"[CELLO] Predicting localization for {len(sequences)} sequences...")

        results = []
        for header, seq in tqdm(sequences, desc="Analyzing localization"):
            protein_id = header.split()[0]

            features = self._calculate_aa_composition(seq)
            scores = self._calculate_scores(features)

            best_loc = max(scores.keys(), key=lambda k: scores[k])
            confidence = scores[best_loc]

            results.append({
                'protein_id': protein_id,
                'description': header,
                'length': len(seq),
                'localization': best_loc,
                'confidence': confidence,
                'cytoplasmic_score': scores.get('Cytoplasmic', 0),
                'inner_membrane_score': scores.get('InnerMembrane', 0),
                'periplasmic_score': scores.get('Periplasmic', 0),
                'outer_membrane_score': scores.get('OuterMembrane', 0),
                'extracellular_score': scores.get('Extracellular', 0),
            })

        return pd.DataFrame(results)

    def run(
        self,
        input_file: str,
        output_dir: str,
    ) -> pd.DataFrame:
        """
        Run localization prediction pipeline.

        Args:
            input_file: Input FASTA file
            output_dir: Output directory

        Returns:
            DataFrame with predictions
        """
        ensure_dir(output_dir)

        sequences = read_fasta(input_file)
        print(f"[CELLO] Loaded {len(sequences)} sequences from {input_file}")

        results = self.predict(sequences)

        output_file = Path(output_dir) / "cello_predictions.tsv"
        results.to_csv(output_file, sep='\t', index=False)
        print(f"[CELLO] Results saved to {output_file}")

        print(f"\n[CELLO] Localization summary:")
        for loc in self.LOCALIZATIONS:
            count = (results['localization'] == loc).sum()
            print(f"  {loc}: {count}")

        return results
