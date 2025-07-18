# CancerPAM: Tumor-Specific Off-Target Analysis Pipeline

## Overview

CancerPAM is a comprehensive pipeline for analyzing tumor-specific off-target risks in CRISPR-Cas9 gene editing. It estimates the probability that tumor-specific mutations could create new off-target sites or rescue existing ones.

## Features

### **Simple Python-Implementation**
   - No external dependencies beyond Python packages
   - Pure Python sequence scanning
   - Fast and reliable
   - Perfect for research and analysis


### 📊 **Analysis Capabilities**

- **Off-target detection**: Find sequences with 0-4 mismatches
- **PAM analysis**: Identify NGG PAM-adjacent sites
- **Rescue risk**: Calculate probability of 4-mismatch → 3-mismatch conversion
- **Novel PAM risk**: Estimate new PAM creation from mutations
- **Comprehensive reporting**: Detailed risk assessment with interpretations

## Quick Start

### 1. Installation

```bash
# Clone or download the scripts
git clone <repository-url>
cd cancerpam

# Install required packages
pip install numpy tqdm pyfaidx biopython

# Run setup (choose option 1 for simple version)
python setup_cancerpam_tools.py
```

### 2. Basic Usage

```bash
# Run analysis with your gRNA and tumor SNV count
python simple_cancerpam.py --gRNA CCGCCAGCGCCGTCTACGTG --tumor_snvs 60000

# Save results to file
python simple_cancerpam.py --gRNA CCGCCAGCGCCGTCTACGTG --tumor_snvs 60000 --output results.txt

# Use custom genome
python simple_cancerpam.py --gRNA CCGCCAGCGCCGTCTACGTG --tumor_snvs 60000 --genome /path/to/genome.fa
```

### 3. Test the Installation

```bash
# Run test suite
python test_cancerpam.py
```

## Detailed Setup Instructions

This is the easiest option and sufficient for most analyses:

1. Install Python 3.7+
2. Install required packages:
   ```bash
   pip install numpy tqdm pyfaidx biopython pandas
   ```
3. Run setup script:
   ```bash
   python setup_cancerpam_tools.py
   ```
   Choose option **1** when prompted.

4. Test installation:
   ```bash
   python test_cancerpam.py
   ```


## Usage Examples

### Basic Analysis
```bash
python simple_cancerpam.py --gRNA CCGCCAGCGCCGTCTACGTG --tumor_snvs 60000
```

### Custom Parameters
```bash
python simple_cancerpam.py \
  --gRNA CCGCCAGCGCCGTCTACGTG \
  --tumor_snvs 60000 \
  --max_mismatches 3 \
  --output detailed_results.txt
```

### Using Full Genome
```bash
# Download full genome first
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz

# Run analysis
python simple_cancerpam.py \
  --gRNA CCGCCAGCGCCGTCTACGTG \
  --tumor_snvs 60000 \
  --genome hg38.fa.gz
```

## Input Parameters

| Parameter | Description | Required | Default |
|-----------|-------------|----------|---------|
| `--gRNA` | 20-nucleotide guide RNA sequence (ACGT only) | Yes | - |
| `--tumor_snvs` | Number of tumor-specific SNVs | Yes | - |
| `--max_mismatches` | Maximum mismatches to analyze (1-4) | No | 4 |
| `--genome` | Path to reference genome FASTA file | No | chr22.fa.gz |
| `--output` | Output file for results | No | stdout |
| `--no_multiprocessing` | Disable parallel processing | No | False |

## Output Interpretation

The pipeline provides comprehensive risk assessment:

### Risk Levels
- 🟢 **LOW** (<1%): Very low risk of tumor-specific off-targets
- 🟡 **MODERATE** (1-5%): Moderate risk - consider additional validation
- 🔴 **HIGH** (>5%): High risk - strong recommendation for additional screening

### Key Metrics
- **PAM-adjacent sites**: Off-targets with NGG PAM by mismatch count
- **Rescue risk**: Probability of 4-mismatch sites becoming active
- **Novel PAM risk**: Probability of new PAM creation near matches
- **Combined risk**: Overall probability of tumor-specific off-targets

## Advanced Features

### Multiprocessing
The pipeline automatically uses multiple CPU cores for faster analysis:
```bash
# Disable multiprocessing (for debugging)
python simple_cancerpam.py --gRNA CCGCCAGCGCCGTCTACGTG --tumor_snvs 60000 --no_multiprocessing
```

### Custom Genome Files
You can use any reference genome in FASTA format:
```bash
python simple_cancerpam.py --gRNA CCGCCAGCGCCGTCTACGTG --tumor_snvs 60000 --genome /path/to/custom_genome.fa
```

## Troubleshooting

### Common Issues

1. **Import errors**: Make sure all required packages are installed
   ```bash
   pip install numpy tqdm pyfaidx biopython
   ```

2. **Memory issues**: Use chromosome-by-chromosome analysis for large genomes
   ```bash
   # The pipeline automatically handles this
   ```

3. **Slow performance**: 
   - Use the sample chromosome (chr22) for testing
   - Enable multiprocessing (default)
   - Consider using the enhanced version with alignment tools

### Getting Help

1. Run the test suite:
   ```bash
   python test_cancerpam.py
   ```

2. Check the configuration:
   ```bash
   ls CancerPAM_tools/
   cat CancerPAM_tools/simple_config.py
   ```

3. Try with the sample data:
   ```bash
   python simple_cancerpam.py --gRNA CCGCCAGCGCCGTCTACGTG --tumor_snvs 60000
   ```

## Algorithm Details

### Off-target Detection
1. **Sequence Scanning**: Exhaustive search for gRNA-like sequences
2. **Mismatch Counting**: Count differences between gRNA and target
3. **PAM Identification**: Check for NGG PAM adjacent to targets
4. **Strand Consideration**: Analyze both forward and reverse strands

### Risk Calculation
1. **Rescue Risk**: Based on Poisson model for SNV rescue events
2. **Novel PAM Risk**: Probability of creating new NGG PAMs
3. **Combined Risk**: Independent probabilities combined

### Performance Optimization
- Multiprocessing for chromosome-level parallelization
- Efficient string matching algorithms
- Memory-efficient genome loading

## File Structure

```
CancerPAM_tools/
├── simple_config.py          # Configuration file
├── chr22.fa.gz               # Sample genome (chromosome 22)

Scripts:
├── simple_cancerpam.py       # Main pipeline (simple version)
├── setup_cancerpam_tools.py  # Setup script
├── test_cancerpam.py         # Test suite
```

## Citation

If you use CancerPAM in your research, please cite:

```
[to be added]
```

## Support

For questions, issues, or contributions:
- Create an issue on GitHub
- Email: [Janik.Hohmann@med.uni-duesseldorf.de]
- Documentation: [link to documentation]
