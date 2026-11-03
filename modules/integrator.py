#!/usr/bin/env python3
"""
Results Integration Module

Combines predictions from all bacterial analysis tools into a unified report.
"""

import os
from pathlib import Path
from typing import List, Dict, Optional
import pandas as pd
import numpy as np

from common import ensure_dir


class ResultsIntegrator:
    """
    Results Integrator

    Combines SignalP, CELLO, and Transmembrane predictions
    to identify surface-exposed proteins.
    """

    SURFACE_LOCALIZATIONS = [
        'OuterMembrane',
        'Extracellular',
        'Periplasmic',
        'CellWall',
    ]

    def __init__(self):
        """Initialize Results Integrator."""
        pass

    def integrate(
        self,
        signalp_df: Optional[pd.DataFrame] = None,
        cello_df: Optional[pd.DataFrame] = None,
        tm_df: Optional[pd.DataFrame] = None,
    ) -> pd.DataFrame:
        """
        Integrate predictions from all tools.

        Args:
            signalp_df: SignalP predictions DataFrame
            cello_df: CELLO predictions DataFrame
            tm_df: Transmembrane predictions DataFrame

        Returns:
            Integrated DataFrame
        """
        print("[Integrator] Combining predictions...")

        dfs = [df for df in [signalp_df, cello_df, tm_df] if df is not None]
        if not dfs:
            raise ValueError("At least one prediction DataFrame is required")

        all_proteins = set()
        for df in dfs:
            all_proteins.update(df['protein_id'].tolist())

        print(f"[Integrator] Total proteins: {len(all_proteins)}")

        results = []
        for protein_id in sorted(all_proteins):
            row = {'protein_id': protein_id}

            # SignalP data
            if signalp_df is not None:
                sp_row = signalp_df[signalp_df['protein_id'] == protein_id]
                if len(sp_row) > 0:
                    sp_row = sp_row.iloc[0]
                    row['has_signal_peptide'] = sp_row.get('has_signal_peptide', False)
                    row['sp_score'] = sp_row.get('sp_score', 0)
                    row['cleavage_site'] = sp_row.get('cleavage_site', 0)
                else:
                    row['has_signal_peptide'] = False
                    row['sp_score'] = 0
                    row['cleavage_site'] = 0

            # CELLO data
            if cello_df is not None:
                cello_row = cello_df[cello_df['protein_id'] == protein_id]
                if len(cello_row) > 0:
                    cello_row = cello_row.iloc[0]
                    row['localization'] = cello_row.get('localization', 'Unknown')
                    row['loc_confidence'] = cello_row.get('confidence', 0)
                    row['length'] = cello_row.get('length', 0)
                else:
                    row['localization'] = 'Unknown'
                    row['loc_confidence'] = 0
                    row['length'] = 0

            # Transmembrane data
            if tm_df is not None:
                tm_row = tm_df[tm_df['protein_id'] == protein_id]
                if len(tm_row) > 0:
                    tm_row = tm_row.iloc[0]
                    row['tm_count'] = tm_row.get('tm_count', 0)
                    row['tm_classification'] = tm_row.get('classification', 'Unknown')
                    row['is_membrane_protein'] = tm_row.get('is_membrane_protein', False)
                else:
                    row['tm_count'] = 0
                    row['tm_classification'] = 'Unknown'
                    row['is_membrane_protein'] = False

            # Determine if surface protein
            is_surface = False
            loc = row.get('localization', '')
            if loc in self.SURFACE_LOCALIZATIONS:
                is_surface = True
            if row.get('has_signal_peptide', False):
                is_surface = True
            if row.get('tm_classification', '') == 'BETA':
                is_surface = True

            row['is_surface_protein'] = is_surface

            # Calculate confidence score
            confidence = 0.0
            if is_surface:
                if loc in self.SURFACE_LOCALIZATIONS:
                    confidence += row.get('loc_confidence', 0) * 0.5
                if row.get('has_signal_peptide', False):
                    confidence += row.get('sp_score', 0) * 0.3
                if row.get('tm_classification', '') == 'BETA':
                    confidence += 0.2
                confidence = min(confidence, 1.0)

            row['surface_confidence'] = confidence
            results.append(row)

        df = pd.DataFrame(results)

        column_order = [
            'protein_id', 'length', 'is_surface_protein', 'surface_confidence',
            'localization', 'loc_confidence', 'has_signal_peptide', 'sp_score',
            'cleavage_site', 'tm_count', 'tm_classification', 'is_membrane_protein',
        ]
        df = df[[c for c in column_order if c in df.columns]]

        return df

    def run(
        self,
        output_dir: str,
        signalp_file: Optional[str] = None,
        cello_file: Optional[str] = None,
        tm_file: Optional[str] = None,
    ) -> pd.DataFrame:
        """
        Run integration pipeline.

        Args:
            output_dir: Output directory
            signalp_file: Path to SignalP results
            cello_file: Path to CELLO results
            tm_file: Path to transmembrane results

        Returns:
            Integrated DataFrame
        """
        ensure_dir(output_dir)

        signalp_df = None
        cello_df = None
        tm_df = None

        if signalp_file and os.path.exists(signalp_file):
            signalp_df = pd.read_csv(signalp_file, sep='\t')
            print(f"[Integrator] Loaded SignalP: {len(signalp_df)} proteins")

        if cello_file and os.path.exists(cello_file):
            cello_df = pd.read_csv(cello_file, sep='\t')
            print(f"[Integrator] Loaded CELLO: {len(cello_df)} proteins")

        if tm_file and os.path.exists(tm_file):
            tm_df = pd.read_csv(tm_file, sep='\t')
            print(f"[Integrator] Loaded TM: {len(tm_df)} proteins")

        results = self.integrate(signalp_df, cello_df, tm_df)

        output_file = Path(output_dir) / "surface_proteins.tsv"
        results.to_csv(output_file, sep='\t', index=False)
        print(f"[Integrator] Results saved to {output_file}")

        surface_count = results['is_surface_protein'].sum()
        print(f"\n[Integrator] Summary:")
        print(f"  Total proteins: {len(results)}")
        print(f"  Surface proteins: {surface_count}")
        print(f"  Percentage: {100*surface_count/len(results):.1f}%")

        high_conf = results[results['surface_confidence'] > 0.5]
        print(f"  High confidence (>0.5): {len(high_conf)}")

        return results
