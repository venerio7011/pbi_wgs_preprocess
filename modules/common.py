#!/usr/bin/env python3
"""
Common Utility Functions

Shared utilities for file handling, environment setup, and device configuration.
"""

import os
from pathlib import Path
from typing import List, Tuple, Optional
import numpy as np


def read_fasta(file_path: str) -> List[Tuple[str, str]]:
    """
    Read sequences from a FASTA file.

    Args:
        file_path: Path to FASTA file

    Returns:
        List of (header, sequence) tuples
    """
    sequences = []
    current_header = None
    current_seq = []

    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith('>'):
                if current_header is not None:
                    sequences.append((current_header, ''.join(current_seq)))
                current_header = line[1:]
                current_seq = []
            else:
                current_seq.append(line)

        if current_header is not None:
            sequences.append((current_header, ''.join(current_seq)))

    return sequences


def write_fasta(
    sequences: List[Tuple[str, str]],
    file_path: str,
    line_width: int = 60
) -> int:
    """
    Write sequences to a FASTA file.

    Args:
        sequences: List of (header, sequence) tuples
        file_path: Path to output file
        line_width: Number of characters per line

    Returns:
        Number of sequences written
    """
    with open(file_path, 'w') as f:
        for header, seq in sequences:
            if not header.startswith('>'):
                header = '>' + header
            f.write(f"{header}\n")
            for i in range(0, len(seq), line_width):
                f.write(f"{seq[i:i+line_width]}\n")

    return len(sequences)


def ensure_dir(directory: str) -> Path:
    """
    Ensure a directory exists, creating it if necessary.
    """
    path = Path(directory)
    path.mkdir(parents=True, exist_ok=True)
    return path


def get_device() -> str:
    """
    Get the best available device for computation.
    """
    try:
        import torch
        if torch.cuda.is_available():
            return 'cuda'
        elif hasattr(torch.backends, 'mps') and torch.backends.mps.is_available():
            return 'mps'
        else:
            return 'cpu'
    except ImportError:
        return 'cpu'


def setup_environment(
    output_dir: Optional[str] = None,
    seed: int = 42,
    verbose: bool = True
) -> dict:
    """
    Setup the analysis environment.
    """
    import random

    random.seed(seed)
    np.random.seed(seed)

    try:
        import torch
        torch.manual_seed(seed)
        if torch.cuda.is_available():
            torch.cuda.manual_seed_all(seed)
    except ImportError:
        pass

    device = get_device()

    if output_dir:
        ensure_dir(output_dir)

    env_info = {
        'device': device,
        'seed': seed,
        'output_dir': output_dir,
    }

    if verbose:
        print(f"[Setup] Device: {device}")
        print(f"[Setup] Random seed: {seed}")
        if output_dir:
            print(f"[Setup] Output directory: {output_dir}")

    return env_info


def get_sequence_stats(sequences: List[Tuple[str, str]]) -> dict:
    """
    Calculate basic statistics for sequences.
    """
    if not sequences:
        return {'count': 0}

    lengths = [len(seq) for _, seq in sequences]

    return {
        'count': len(sequences),
        'total_length': sum(lengths),
        'min_length': min(lengths),
        'max_length': max(lengths),
        'mean_length': np.mean(lengths),
        'median_length': np.median(lengths),
    }
