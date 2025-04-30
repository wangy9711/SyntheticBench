# Benchmarking Full-Atom Structure Prediction for Peptide–Protein Complexes

## Description

This repository contains the evaluation code for our benchmark study on full-atom peptide–protein structure prediction algorithms. We systematically evaluated multiple state-of-the-art models (e.g., AlphaFold 2, AlphaFold 3, ProteniX, Chai-1, HelixFold 3) on a curated dataset of canonical, non-canonical, and cyclic peptides.

We assessed model performance across different input modalities—including amino acid sequences, SMILES strings, post-translational modification (PTM) annotations, and covalent bond constraints—using quantitative metrics such as RMSD, center-of-mass deviation, and binding site clustering success rates. Our evaluation covers both the global structure of the complex and the intrinsic conformational accuracy of the peptide component.

## Workflow Overview

![Workflow Diagram](./workflow.png)

## Repository Structure

```
/analysis
  ├── add_clust.py     # Add the number of conformational clusters after final scoring
  ├── align_save.py    # Output conformations after alignment 
  ├── cal_score.py     # Calculate all scores for prediction
  └── utils.py         # Toolkit


/visualization
  └── draw.ipynb       # Draw pictures based on the scores

```

## Dependencies

- Python ≥ 3.8
- numpy
- scipy
- rdkit
- scikit-learn
- seaborn
- matplotlib
- pymol
- networkx
- biopython

## Data
All datasets used in this work will be available once the manuscript is accepted for publication.



