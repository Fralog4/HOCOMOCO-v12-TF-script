# DNA Binding Site Analyzer

This Python program is designed to analyze DNA sequences and identify binding sites for a given Transcription Factor (TF) using a Position Weight Matrix (PWM).

## Features

- Loads a DNA sequence from a FASTA file
- Loads a TF binding site PWM from a Transfac format file
- Scans the DNA sequence to find binding sites that match the PWM
- Visualizes the distribution of binding site scores using a histogram
- Prints the position and score of each binding site found

## Requirements

This program requires the following Python libraries:

- Biopython
- Matplotlib
- Numpy

You can install the required dependencies using pip:

```bash
pip install biopython matplotlib
``` 