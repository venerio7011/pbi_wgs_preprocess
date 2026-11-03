#!/usr/bin/env python3
"""
Prediction Merger Module

Merge CELLO and SignalP 6.0 predictions into a unified output.
"""

import sys
from pathlib import Path
from typing import Dict, List, Any


def parse_cello_results(file_path: str) -> Dict[str, Dict[str, Any]]:
    """
    Parse CELLO output file.

    Args:
        file_path: Path to CELLO results file

    Returns:
        Dictionary mapping protein IDs to predictions
    """
    predictions = {}

    if not Path(file_path).exists():
        print(f"Warning: CELLO results not found at {file_path}")
        return predictions

    with open(file_path, 'r') as f:
        in_results = False
        for line in f:
            line = line.strip()

            if line.startswith('SeqID'):
                in_results = True
                continue

            if in_results and line and not line.startswith('-') and not line.startswith('='):
                if line.startswith('Summary:'):
                    break

                parts = line.split()
                if len(parts) >= 3:
                    seq_id = parts[0]
                    location = parts[2]
                    score = parts[3] if len(parts) > 3 else "N/A"
                    predictions[seq_id] = {'location': location, 'score': score}

    return predictions


def parse_signalp6_results(file_path: str) -> Dict[str, Dict[str, Any]]:
    """
    Parse SignalP 6.0 output file with all details.

    Args:
        file_path: Path to SignalP 6.0 results file

    Returns:
        Dictionary mapping protein IDs to predictions
    """
    predictions = {}

    if not Path(file_path).exists():
        print(f"Warning: SignalP 6.0 results not found at {file_path}")
        return predictions

    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            parts = line.split('\t')
            if len(parts) >= 2:
                full_id = parts[0]
                seq_id = full_id.split()[0] if ' ' in full_id else full_id

                prediction = parts[1] if len(parts) > 1 else "OTHER"
                other_score = parts[2] if len(parts) > 2 else "N/A"
                sp_score = parts[3] if len(parts) > 3 else "N/A"
                lipo_score = parts[4] if len(parts) > 4 else "N/A"
                tat_score = parts[5] if len(parts) > 5 else "N/A"
                tatlipo_score = parts[6] if len(parts) > 6 else "N/A"
                pilin_score = parts[7] if len(parts) > 7 else "N/A"
                cs_position = parts[8].strip() if len(parts) > 8 else "N/A"

                predictions[seq_id] = {
                    'prediction': prediction,
                    'other_score': other_score,
                    'sp_score': sp_score,
                    'lipo_score': lipo_score,
                    'tat_score': tat_score,
                    'tatlipo_score': tatlipo_score,
                    'pilin_score': pilin_score,
                    'cs_position': cs_position
                }

    return predictions


def merge_predictions(
    cello_file: str,
    signalp6_file: str,
    output_file: str
) -> List[Dict[str, Any]]:
    """
    Merge predictions from CELLO and SignalP 6.0.

    Args:
        cello_file: Path to CELLO results
        signalp6_file: Path to SignalP 6.0 results
        output_file: Path to output file

    Returns:
        List of merged prediction dictionaries
    """
    print("Parsing CELLO results...")
    cello_pred = parse_cello_results(cello_file)

    print("Parsing SignalP 6.0 results...")
    signalp6_pred = parse_signalp6_results(signalp6_file)

    all_ids = set(cello_pred.keys()) | set(signalp6_pred.keys())

    print(f"\nTotal unique proteins: {len(all_ids)}")
    print(f"CELLO predictions: {len(cello_pred)}")
    print(f"SignalP 6.0 predictions: {len(signalp6_pred)}")

    merged = []

    for seq_id in sorted(all_ids):
        cello_info = cello_pred.get(seq_id, {'location': 'N/A', 'score': 'N/A'})
        cello_loc = cello_info['location']
        cello_score = cello_info.get('score', 'N/A')

        signalp6_info = signalp6_pred.get(seq_id, {
            'prediction': 'N/A',
            'other_score': 'N/A',
            'sp_score': 'N/A',
            'lipo_score': 'N/A',
            'tat_score': 'N/A',
            'tatlipo_score': 'N/A',
            'pilin_score': 'N/A',
            'cs_position': 'N/A'
        })

        merged.append({
            'id': seq_id,
            'sp_pred': signalp6_info['prediction'],
            'sp_other': signalp6_info['other_score'],
            'sp_sp': signalp6_info['sp_score'],
            'sp_lipo': signalp6_info['lipo_score'],
            'sp_tat': signalp6_info['tat_score'],
            'sp_tatlipo': signalp6_info['tatlipo_score'],
            'sp_pilin': signalp6_info['pilin_score'],
            'sp_cs': signalp6_info['cs_position'],
            'cello': cello_loc,
            'cello_score': cello_score
        })

    with open(output_file, 'w') as f:
        f.write("=" * 100 + "\n")
        f.write("Merged Protein Subcellular Localization Predictions\n")
        f.write("CELLO + SignalP 6.0 (Streamlined Pipeline)\n")
        f.write("=" * 100 + "\n\n")

        f.write(f"Total proteins analyzed: {len(all_ids)}\n")
        f.write(f"SignalP 6.0 predictions: {len(signalp6_pred)}\n")
        f.write(f"CELLO predictions: {len(cello_pred)}\n")

        proteins_with_sp = sum(1 for m in merged if m.get('sp_pred') not in ['OTHER', 'N/A'])
        f.write(f"Proteins with signal peptides: {proteins_with_sp}\n")

        f.write("\n" + "=" * 140 + "\n\n")

        header = "{:<15} {:<10} {:<10} {:<10} {:<10} {:<10} {:<10} {:<10} {:<10} {:<18} {:<8}\n"
        f.write(header.format(
            'Protein ID', 'SP_Pred', 'OTHER', 'SP', 'LIPO', 'TAT', 'TATLIPO', 'PILIN', 'CS_Pos', 'CELLO', 'CELLO_Conf'))
        f.write("-" * 140 + "\n")

        for entry in merged:
            f.write(header.format(
                entry['id'],
                entry['sp_pred'],
                entry['sp_other'],
                entry['sp_sp'],
                entry['sp_lipo'],
                entry['sp_tat'],
                entry['sp_tatlipo'],
                entry['sp_pilin'],
                entry['sp_cs'],
                entry['cello'],
                entry['cello_score']))

    print(f"\nMerged results written to: {output_file}")
    return merged


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: prediction_merger.py <cello_file> <signalp6_file> <output_file>")
        sys.exit(1)

    cello_file = sys.argv[1]
    signalp6_file = sys.argv[2]
    output_file = sys.argv[3]

    merge_predictions(cello_file, signalp6_file, output_file)
