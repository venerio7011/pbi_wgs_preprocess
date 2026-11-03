# PBI Unified Pipeline

**Phage-Bacteria Interaction (PBI) Analysis Pipeline**

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/venerio7011/pbi_wgs_preprocess/blob/main/PBI_Unified_Pipeline.ipynb)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)

A comprehensive Google Colab-ready pipeline for analyzing phage-bacteria interactions by identifying:
1. **Bacterial surface proteins** (potential phage receptors)
2. **Phage receptor-binding proteins (RBPs)** and **depolymerases**

## Features

### Part A: Bacterial Surface Protein Analysis
| Module | Function | Method |
|--------|----------|--------|
| SignalP-like | Signal peptide prediction | N-terminal sequence features |
| CELLO-like | Subcellular localization | Amino acid composition |
| Transmembrane | TM topology prediction | Hydrophobicity profile |
| Integration | Surface protein classification | Multi-evidence scoring |

### Part B: Phage RBP/DPO Analysis
| Module | Function | Method |
|--------|----------|--------|
| RBP Detection | Receptor-binding protein identification | Sequence features + keywords |
| DPO Prediction | Depolymerase detection | Catalytic residue analysis |
| Integration | Combined classification | Multi-score integration |

## Quick Start

### Option 1: Google Colab (Recommended)

1. Click the "Open in Colab" badge above
2. Upload your bacterial and phage protein FASTA files
3. Run all cells sequentially

### Option 2: Local Installation

```bash
# Clone repository
git clone https://github.com/venerio7011/pbi_wgs_preprocess.git
cd pbi_wgs_preprocess

# Create virtual environment (recommended)
python -m venv venv
source venv/bin/activate  # Linux/macOS
# or: venv\Scripts\activate  # Windows

# Install dependencies
pip install -r requirements.txt

# Run Jupyter notebook
jupyter notebook PBI_Unified_Pipeline.ipynb
```

## Input Requirements

| File | Description | Format |
|------|-------------|--------|
| Bacterial proteins | Proteome FASTA | `.faa` or `.fasta` |
| Phage proteins | Proteome FASTA | `.faa` or `.fasta` |

### Input File Format

```
>protein_id description
MKFLILLFNILCLFPVLAADNHGVGPGACGNLLGGLTDPFDGLAGRKLVMVGDGLDVYTP...
>protein_id_2 another description
MKTAIAIAVALAGFATVAQAAPKDNTWYTGAKLGWSQYHDTGFINNNGPTHENQLGAGA...
```

## Output Structure

```
output/
├── bacteria/
│   ├── signalp/
│   │   └── signalp_predictions.tsv
│   ├── cello/
│   │   └── cello_predictions.tsv
│   ├── transmembrane/
│   │   └── transmembrane_predictions.tsv
│   └── integrated/
│       └── surface_proteins.tsv
├── phage/
│   ├── rbpdetect2/
│   │   └── rbp_predictions.tsv
│   ├── phagedpo/
│   │   └── depolymerase_predictions.tsv
│   └── integrated/
│       └── rbp_dpo_predictions.tsv
└── pbi_analysis/
    ├── high_conf_surface_proteins.tsv
    ├── high_conf_rbp_dpo.tsv
    └── pbi_summary.png
```

## Pipeline Overview

```
┌─────────────────────────────────────────────────────────────┐
│                    PBI Unified Pipeline                      │
├─────────────────────────────┬───────────────────────────────┤
│   BACTERIAL ANALYSIS        │      PHAGE ANALYSIS           │
├─────────────────────────────┼───────────────────────────────┤
│                             │                               │
│  ┌─────────────────┐        │   ┌─────────────────┐         │
│  │    SignalP      │        │   │   RBP Detector  │         │
│  │ (Signal peptide)│        │   │ (Tail fibers)   │         │
│  └────────┬────────┘        │   └────────┬────────┘         │
│           │                 │            │                  │
│  ┌────────▼────────┐        │   ┌────────▼────────┐         │
│  │     CELLO       │        │   │  Depolymerase   │         │
│  │ (Localization)  │        │   │   Predictor     │         │
│  └────────┬────────┘        │   └────────┬────────┘         │
│           │                 │            │                  │
│  ┌────────▼────────┐        │            │                  │
│  │  Transmembrane  │        │            │                  │
│  │   Predictor     │        │            │                  │
│  └────────┬────────┘        │            │                  │
│           │                 │            │                  │
│  ┌────────▼────────┐        │   ┌────────▼────────┐         │
│  │   Integration   │        │   │   Integration   │         │
│  └────────┬────────┘        │   └────────┬────────┘         │
│           │                 │            │                  │
├───────────┴─────────────────┴────────────┴──────────────────┤
│                                                             │
│                  Combined PBI Analysis Report               │
│           (Surface proteins ↔ RBPs/Depolymerases)           │
│                                                             │
└─────────────────────────────────────────────────────────────┘
```

## Output Columns

### Bacterial Surface Proteins (`surface_proteins.tsv`)

| Column | Description |
|--------|-------------|
| protein_id | Protein identifier |
| length | Sequence length |
| is_surface_protein | Boolean classification |
| surface_confidence | Confidence score (0-1) |
| localization | Predicted localization |
| has_signal_peptide | Signal peptide detected |
| tm_count | Transmembrane helix count |

### Phage RBP/DPO (`rbp_dpo_predictions.tsv`)

| Column | Description |
|--------|-------------|
| protein_id | Protein identifier |
| length | Sequence length |
| rbp_score | RBP prediction score |
| is_rbp | RBP classification |
| dpo_score | Depolymerase score |
| is_depolymerase | DPO classification |
| classification | Combined class (RBP/DPO/RBP+DPO/Other) |

## Dependencies

```
biopython>=1.80
pandas>=1.5.0
numpy>=1.23.0
scipy>=1.9.0
scikit-learn>=1.1.0
matplotlib>=3.6.0
seaborn>=0.12.0
tqdm>=4.64.0
```

## Limitations

- The signal peptide and localization predictors are simplified implementations based on sequence features
- For publication-quality predictions, consider using the original tools:
  - [SignalP 6.0](https://services.healthtech.dtu.dk/services/SignalP-6.0/)
  - [PSORTb](https://www.psort.org/psortb/)
  - [TMHMM](https://services.healthtech.dtu.dk/services/TMHMM-2.0/)
  - [DeepTMHMM](https://dtu.biolib.com/DeepTMHMM)

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

1. Fork the repository
2. Create your feature branch (`git checkout -b feature/AmazingFeature`)
3. Commit your changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- Inspired by established bioinformatics tools for protein analysis
- Built for phage-bacteria interaction research

## Contact

For questions or issues, please open a GitHub issue.
