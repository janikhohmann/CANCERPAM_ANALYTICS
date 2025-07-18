#!/usr/bin/env python3
"""
Simplified CancerPAM Pipeline - Python Only
==========================================

This script provides a pure Python implementation for tumor-specific
off-target analysis without requiring BWA, Bowtie2, or Cas-OFFinder.

Features:
- Pure Python implementation (no external alignment tools needed)
- Efficient sequence scanning with multiprocessing
- Comprehensive off-target analysis
- Detailed risk assessment
"""

import os
import sys
import argparse
import re
import math
import gzip
import urllib.request
from pathlib import Path
from collections import defaultdict
from typing import Dict, List, Tuple
from concurrent.futures import ProcessPoolExecutor, as_completed
import time

# Check for required packages
try:
    import numpy as np
    from tqdm import tqdm
except ImportError as e:
    print(f"❌ Missing package: {e}")
    print("Please install: pip install numpy tqdm")
    sys.exit(1)

class SimplePAMAnalyzer:
    """
    Simple, efficient off-target analyzer using Python only.
    """
    
    def __init__(self, genome_fasta: str):
        self.genome_fasta = Path(genome_fasta)
        self.chromosomes = {}
        self.load_genome()
    
    def load_genome(self):
        """Load reference genome from FASTA file."""
        print("📖 Loading reference genome...")
        print(f"   Genome file: {self.genome_fasta}")
        
        if not self.genome_fasta.exists():
            print(f"❌ Genome file not found: {self.genome_fasta}")
            
            # Ask user what to download
            print("\n🤔 What would you like to download?")
            print("1. Complete human genome (hg38) - ~900 MB download, all chromosomes")
            print("2. Sample chromosome (chr22) - ~12 MB download, for testing")
            
            choice = input("\nChoose option (1 or 2) [1]: ").strip() or "1"
            
            if choice == "1":
                print("Downloading complete human genome...")
                self.download_full_genome()
            else:
                print("Downloading sample chromosome for testing...")
                self.download_sample_genome()
        
        # Check file size
        file_size = self.genome_fasta.stat().st_size
        print(f"   File size: {file_size:,} bytes ({file_size/1024/1024:.1f} MB)")
        
        # Load FASTA file
        current_chrom = None
        current_seq = []
        lines_processed = 0
        
        opener = gzip.open if str(self.genome_fasta).endswith('.gz') else open
        mode = 'rt' if str(self.genome_fasta).endswith('.gz') else 'r'
        
        print(f"   File format: {'gzip compressed' if str(self.genome_fasta).endswith('.gz') else 'plain text'}")
        print("   Processing FASTA file...")
        
        with opener(self.genome_fasta, mode) as f:
            for line in tqdm(f, desc="Loading genome"):
                line = line.strip()
                lines_processed += 1
                
                if line.startswith('>'):
                    # Save previous chromosome
                    if current_chrom and current_seq:
                        seq_length = len(''.join(current_seq))
                        print(f"   ✅ Loaded {current_chrom}: {seq_length:,} bp")
                        self.chromosomes[current_chrom] = ''.join(current_seq).upper()
                    
                    # Start new chromosome
                    current_chrom = line[1:].split()[0]  # Take first part of header
                    print(f"   📍 Starting chromosome: {current_chrom}")
                    current_seq = []
                elif current_chrom:
                    current_seq.append(line)
                
                # Progress report every 10000 lines
                if lines_processed % 10000 == 0:
                    current_size = len(''.join(current_seq)) if current_seq else 0
                    print(f"      Lines processed: {lines_processed:,}, Current chrom size: {current_size:,} bp")
            
            # Save last chromosome
            if current_chrom and current_seq:
                seq_length = len(''.join(current_seq))
                print(f"   ✅ Loaded {current_chrom}: {seq_length:,} bp")
                self.chromosomes[current_chrom] = ''.join(current_seq).upper()
        
        print(f"✅ Genome loading complete!")
        print(f"   Total lines processed: {lines_processed:,}")
        print(f"   Chromosomes loaded: {len(self.chromosomes)}")
        
        total_bp = 0
        for chrom, seq in self.chromosomes.items():
            chrom_size = len(seq)
            total_bp += chrom_size
            print(f"   📍 {chrom}: {chrom_size:,} bp")
        
        print(f"   Total genome size: {total_bp:,} bp ({total_bp/1e6:.1f} Mbp)")
        
        # Validate sequences
        print("   Validating sequences...")
        for chrom, seq in self.chromosomes.items():
            n_count = seq.count('N')
            gc_count = seq.count('G') + seq.count('C')
            at_count = seq.count('A') + seq.count('T')
            gc_content = (gc_count / len(seq)) * 100 if len(seq) > 0 else 0
            
            print(f"   📊 {chrom}: GC content = {gc_content:.1f}%, N's = {n_count:,} ({n_count/len(seq)*100:.2f}%)")
            
            # Check for unusual characters
            valid_bases = set('ACGTN')
            seq_bases = set(seq)
            invalid_bases = seq_bases - valid_bases
            if invalid_bases:
                print(f"   ⚠️  {chrom}: Found unusual characters: {invalid_bases}")
        
        print("   ✅ Sequence validation complete!")
    
    def download_sample_genome(self):
        """Download a sample chromosome for testing."""
        sample_url = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr22.fa.gz"
        sample_path = self.genome_fasta.parent / "chr22.fa.gz"
        
        print(f"⬇️  Downloading sample genome (chr22) to {sample_path}")
        try:
            urllib.request.urlretrieve(sample_url, sample_path)
            self.genome_fasta = sample_path
            print("✅ Sample genome downloaded")
        except Exception as e:
            print(f"❌ Failed to download sample genome: {e}")
            sys.exit(1)
    
    def download_full_genome(self):
        """Download the complete human genome."""
        full_genome_url = "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"
        full_genome_path = self.genome_fasta.parent / "hg38.fa.gz"
        
        print(f"⬇️  Downloading complete human genome (hg38) to {full_genome_path}")
        print("   ⚠️  This is a large file (~900 MB compressed, ~3.1 GB uncompressed)")
        print("   ⏱️  Download may take several minutes depending on your internet connection")
        
        try:
            # Download with progress bar
            def show_progress(block_num, block_size, total_size):
                if total_size > 0:
                    percent = min(100, (block_num * block_size * 100) / total_size)
                    downloaded = block_num * block_size
                    print(f"\r   Progress: {percent:.1f}% ({downloaded:,} / {total_size:,} bytes)", end='')
            
            urllib.request.urlretrieve(full_genome_url, full_genome_path, reporthook=show_progress)
            print("\n✅ Complete genome downloaded")
            self.genome_fasta = full_genome_path
            
        except Exception as e:
            print(f"\n❌ Failed to download complete genome: {e}")
            print("   Falling back to sample chromosome (chr22)...")
            self.download_sample_genome()
    
    
    def reverse_complement(self, seq: str) -> str:
        """Get reverse complement of sequence."""
        complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
        return ''.join(complement.get(base, base) for base in reversed(seq))
    
    def count_mismatches(self, seq1: str, seq2: str) -> int:
        """Count mismatches between two sequences."""
        if len(seq1) != len(seq2):
            return len(seq1)  # Return max if lengths differ
        return sum(1 for a, b in zip(seq1, seq2) if a != b)
    
    def has_ngg_pam(self, chrom_seq: str, pos: int, strand: str, grna_len: int) -> bool:
        """Check if position has NGG PAM."""
        try:
            if strand == '+':
                # PAM is immediately after the gRNA
                pam_start = pos + grna_len
                if pam_start + 2 < len(chrom_seq):
                    pam = chrom_seq[pam_start:pam_start + 3]
                    return pam.endswith('GG')
            else:
                # PAM is before the gRNA on reference (for reverse strand)
                pam_start = pos - 3
                if pam_start >= 0:
                    pam = chrom_seq[pam_start:pam_start + 3]
                    pam_rc = self.reverse_complement(pam)
                    return pam_rc.endswith('GG')
        except IndexError:
            pass
        return False
    
    def scan_chromosome(self, chrom_name: str, chrom_seq: str, grna: str, 
                       grna_rc: str, max_mismatches: int) -> Tuple[Dict[int, int], List[Tuple]]:
        """Scan a single chromosome for off-targets."""
        pam_counts = defaultdict(int)
        no_pam_sites = []
        grna_len = len(grna)
        
        print(f"      🔍 Scanning {chrom_name} ({len(chrom_seq):,} bp)")
        print(f"         Looking for sequences similar to: {grna}")
        print(f"         Reverse complement: {grna_rc}")
        
        # Progress tracking
        positions_scanned = 0
        forward_matches = 0
        reverse_matches = 0
        
        # Adaptive reporting interval based on chromosome size
        if len(chrom_seq) > 100_000_000:  # Large chromosomes (>100M bp)
            report_interval = len(chrom_seq) // 50  # Report every 2%
        elif len(chrom_seq) > 50_000_000:  # Medium chromosomes (>50M bp)
            report_interval = len(chrom_seq) // 25  # Report every 4%
        else:  # Smaller chromosomes
            report_interval = max(len(chrom_seq) // 10, 1_000_000)  # Report every 10% or 1M bp
        
        # Scan forward strand
        print(f"         Scanning forward strand...")
        for i in range(len(chrom_seq) - grna_len + 1):
            target_seq = chrom_seq[i:i + grna_len]
            mismatches = self.count_mismatches(grna, target_seq)
            
            if mismatches <= max_mismatches:
                forward_matches += 1
                if self.has_ngg_pam(chrom_seq, i, '+', grna_len):
                    pam_counts[mismatches] += 1
                    if mismatches <= 2:  # Only show exact and near-exact matches
                        pam_context = chrom_seq[i:i + grna_len + 3] if i + grna_len + 3 <= len(chrom_seq) else chrom_seq[i:]
                        print(f"           Found PAM site at pos {i:,}: {target_seq} -> {pam_context[-3:]} ({mismatches} mismatches)")
                elif mismatches <= 3:  # Only consider ≤3 mismatches for no-PAM sites
                    no_pam_sites.append((chrom_name, i, '+', mismatches))
                    if mismatches <= 1:  # Only show very close matches
                        print(f"           Found no-PAM site at pos {i:,}: {target_seq} ({mismatches} mismatches)")
            
            positions_scanned += 1
            
            # Progress report
            if positions_scanned % report_interval == 0:
                progress = (positions_scanned / (len(chrom_seq) - grna_len + 1)) * 50  # 50% for forward strand
                #print(f"         Forward progress: {progress:.1f}% ({positions_scanned:,}/{len(chrom_seq) - grna_len + 1:,})")
        
        print(f"         Forward strand complete: {forward_matches} matches found")
        
        # Scan reverse strand
        print(f"         Scanning reverse strand...")
        positions_scanned = 0
        for i in range(len(chrom_seq) - grna_len + 1):
            target_seq = chrom_seq[i:i + grna_len]
            mismatches = self.count_mismatches(grna_rc, target_seq)
            
            if mismatches <= max_mismatches:
                reverse_matches += 1
                if self.has_ngg_pam(chrom_seq, i, '-', grna_len):
                    pam_counts[mismatches] += 1
                    if mismatches <= 2:  # Only show exact and near-exact matches
                        pam_context = chrom_seq[max(0, i-3):i + grna_len] if i >= 3 else chrom_seq[:i + grna_len]
                        print(f"           Found PAM site at pos {i:,}: {target_seq} -> {pam_context[:3]} ({mismatches} mismatches)")
                elif mismatches <= 3:  # Only consider ≤3 mismatches for no-PAM sites
                    no_pam_sites.append((chrom_name, i, '-', mismatches))
                    if mismatches <= 1:  # Only show very close matches
                        print(f"           Found no-PAM site at pos {i:,}: {target_seq} ({mismatches} mismatches)")
            
            positions_scanned += 1
            
            # Progress report
            if positions_scanned % report_interval == 0:
                progress = 50 + (positions_scanned / (len(chrom_seq) - grna_len + 1)) * 50  # 50-100% for reverse strand
                #print(f"         Reverse progress: {progress:.1f}% ({positions_scanned:,}/{len(chrom_seq) - grna_len + 1:,})")
        
        print(f"         Reverse strand complete: {reverse_matches} matches found")
        
        # Summary for this chromosome
        total_pam_sites = sum(pam_counts.values())
        total_no_pam = len(no_pam_sites)
        print(f"         {chrom_name} summary:")
        print(f"           Total matches: {forward_matches + reverse_matches}")
        print(f"           PAM-adjacent sites: {total_pam_sites}")
        print(f"           No-PAM sites: {total_no_pam}")
        
        return dict(pam_counts), no_pam_sites
    
    def analyze_off_targets(self, grna: str, max_mismatches: int = 4, 
                           use_multiprocessing: bool = True) -> Tuple[Dict[int, int], List[Tuple]]:
        """Analyze off-targets across the genome."""
        print(f"🔍 Analyzing off-targets for gRNA: {grna}")
        print(f"   Max mismatches: {max_mismatches}")
        print(f"   Multiprocessing: {'enabled' if use_multiprocessing else 'disabled'}")
        print(f"   Chromosomes to scan: {len(self.chromosomes)}")
        
        # Show chromosome sizes
        total_bp = 0
        for chrom, seq in self.chromosomes.items():
            chrom_size = len(seq)
            total_bp += chrom_size
            print(f"   📍 {chrom}: {chrom_size:,} bp")
        
        print(f"   Total genome size: {total_bp:,} bp")
        
        grna_rc = self.reverse_complement(grna)
        print(f"   gRNA reverse complement: {grna_rc}")
        
        total_pam_counts = defaultdict(int)
        all_no_pam_sites = []
        
        if use_multiprocessing and len(self.chromosomes) > 1:
            print("🔄 Using multiprocessing for chromosome scanning...")
            
            # Use multiprocessing for multiple chromosomes
            with ProcessPoolExecutor() as executor:
                print(f"   Submitting {len(self.chromosomes)} jobs to process pool...")
                futures = {
                    executor.submit(self.scan_chromosome, chrom, seq, grna, grna_rc, max_mismatches): chrom
                    for chrom, seq in self.chromosomes.items()
                }
                
                completed_count = 0
                for future in tqdm(as_completed(futures), total=len(futures), desc="Scanning chromosomes"):
                    chrom = futures[future]
                    completed_count += 1
                    
                    try:
                        print(f"   ✅ Processing results from {chrom} ({completed_count}/{len(futures)})")
                        pam_counts, no_pam_sites = future.result()
                        
                        # Show intermediate results
                        pam_sites_found = sum(pam_counts.values())
                        print(f"      Found {pam_sites_found} PAM-adjacent sites, {len(no_pam_sites)} no-PAM sites")
                        
                        # Merge results
                        for mismatches, count in pam_counts.items():
                            total_pam_counts[mismatches] += count
                        all_no_pam_sites.extend(no_pam_sites)
                        
                        # Show running totals
                        running_total_pam = sum(total_pam_counts.values())
                        running_total_no_pam = len(all_no_pam_sites)
                        print(f"      Running totals: {running_total_pam} PAM sites, {running_total_no_pam} no-PAM sites")
                        
                    except Exception as e:
                        print(f"❌ Error processing {chrom}: {e}")
        else:
            print("🔄 Using sequential processing...")
            
            # Sequential processing
            chrom_count = 0
            for chrom, seq in tqdm(self.chromosomes.items(), desc="Scanning chromosomes"):
                chrom_count += 1
                print(f"   🔍 Scanning chromosome {chrom} ({chrom_count}/{len(self.chromosomes)})")
                print(f"      Length: {len(seq):,} bp")
                
                pam_counts, no_pam_sites = self.scan_chromosome(chrom, seq, grna, grna_rc, max_mismatches)
                
                # Show results for this chromosome
                pam_sites_found = sum(pam_counts.values())
                print(f"      Results: {pam_sites_found} PAM-adjacent sites, {len(no_pam_sites)} no-PAM sites")
                
                if pam_counts:
                    for mismatches, count in sorted(pam_counts.items()):
                        print(f"        {mismatches} mismatches: {count} sites")
                
                # Merge results
                for mismatches, count in pam_counts.items():
                    total_pam_counts[mismatches] += count
                all_no_pam_sites.extend(no_pam_sites)
                
                # Show running totals
                running_total_pam = sum(total_pam_counts.values())
                running_total_no_pam = len(all_no_pam_sites)
                print(f"      Running totals: {running_total_pam} PAM sites, {running_total_no_pam} no-PAM sites")
        
        # Final summary
        print(f"📊 Final results summary:")
        print(f"   Total PAM-adjacent sites: {sum(total_pam_counts.values())}")
        print(f"   Total no-PAM sites: {len(all_no_pam_sites)}")
        print(f"   Breakdown by mismatches:")
        for mismatches in sorted(total_pam_counts.keys()):
            print(f"     {mismatches} mismatches: {total_pam_counts[mismatches]} sites")
        
        return dict(total_pam_counts), all_no_pam_sites

def calculate_rescue_risk(n_offtargets_4: int, n_mutations: int, 
                         grna_length: int = 20, genome_size: float = 3.1e9) -> Tuple[float, float]:
    """Calculate risk of tumor SNVs rescuing 4-mismatch sites to 3-mismatch."""
    # Probability that an SNV corrects a mismatch
    p_correct_mismatch = 1/3 # Assuming 1/3 chance of correcting a mismatch because of random mutation
    
    # Expected number of rescue events
    lambda_rescue = n_offtargets_4 * n_mutations * (4 / genome_size) * p_correct_mismatch
    
    # Probability of at least one rescue event
    prob_rescue = 1 - math.exp(-lambda_rescue)
    
    return lambda_rescue, prob_rescue

def calculate_novel_pam_risk(no_pam_count: int, n_mutations: int,
                           genome_size: float = 3.1e9) -> Tuple[float, float]:
    """Calculate risk of tumor SNVs creating novel NGG PAMs."""
    # Probability that an SNV in PAM region creates NGG
    # 6 out of 15 dinucleotides can mutate to GG (GG excluded beacause it can#t mutate to GG anymore): AG,TG,CG,GA,GT,GC
    # P(dinucleotide can mutate to GG) = 6/15 = 0.4
    # P(mutation creates GG | dinucleotide is mutable) = 1/3
    # P(single strand) = 6/15 * 1/3 = 0.13333333
    p_create_pam = 0.133333  # Probability of creating a novel PAM
    
    # Expected number of novel PAM events
    lambda_pam = no_pam_count * n_mutations * (3 / genome_size) * p_create_pam
    
    # Probability of at least one novel PAM event
    prob_pam = 1 - math.exp(-lambda_pam)
    
    return lambda_pam, prob_pam

def validate_grna(grna: str) -> bool:
    """Validate gRNA sequence."""
    if len(grna) != 20:
        print(f"❌ gRNA must be 20 nucleotides, got {len(grna)}")
        return False
    
    if not re.match(r'^[ACGT]+$', grna.upper()):
        print("❌ gRNA must contain only A, C, G, T nucleotides")
        return False
    
    return True

def generate_detailed_report(grna: str, tumor_snvs: int, pam_counts: Dict[int, int],
                           no_pam_count: int, rescue_lambda: float, rescue_prob: float,
                           novel_pam_lambda: float, novel_pam_prob: float) -> str:
    """Generate comprehensive analysis report."""
    
    total_pam_sites = sum(pam_counts.values())
    
    report = f"""
╔═════════════════════════════════════════════════════════════════════════════════════╗
║                          🧬 CancerPAM Analysis Report                               ║
╚═════════════════════════════════════════════════════════════════════════════════════╝

📊 INPUT PARAMETERS
───────────────────────────────────────────────────────────────────────────────────────
  gRNA sequence:        {grna}
  Tumor SNVs:           {tumor_snvs:,}
  Analysis method:      Pure Python implementation

📈 OFF-TARGET ANALYSIS RESULTS
───────────────────────────────────────────────────────────────────────────────────────
  Sites with NGG PAM:
"""
    
    for mismatches in sorted(pam_counts.keys()):
        count = pam_counts[mismatches]
        percentage = (count / total_pam_sites * 100) if total_pam_sites > 0 else 0
        report += f"    {mismatches} mismatches: {count:,} sites ({percentage:.1f}%)\n"
    
    report += f"""
  Total PAM-adjacent sites: {total_pam_sites:,}
  Sites without NGG PAM:    {no_pam_count:,}

⚠️  RISK ASSESSMENT
───────────────────────────────────────────────────────────────────────────────────────
  1. Rescue Risk (4-mismatch → ≤3-mismatch):
     4-mismatch sites with PAM: {pam_counts.get(4, 0):,}
     Expected rescue events:    {rescue_lambda:.6f}
     Probability:               {rescue_prob:.2%}
     
  2. Novel PAM Creation Risk:
     Near-matches without PAM:  {no_pam_count:,}
     Expected PAM events:       {novel_pam_lambda:.6f}
     Probability:               {novel_pam_prob:.2%}

🎯 COMBINED RISK SUMMARY
───────────────────────────────────────────────────────────────────────────────────────
  Individual risks are independent, so combined probability is:
  P(rescue OR novel PAM) ≈ {(rescue_prob + novel_pam_prob):.2%}
  
  Risk interpretation:
"""
    
    combined_risk = rescue_prob + novel_pam_prob
    if combined_risk < 0.01:
        risk_level = "🟢 LOW"
        interpretation = "Very low risk of tumor-specific off-targets"
    elif combined_risk < 0.05:
        risk_level = "🟡 MODERATE"
        interpretation = "Moderate risk - consider additional validation"
    else:
        risk_level = "🔴 HIGH"
        interpretation = "High risk - strong recommendation for additional screening"
    
    report += f"  Risk level: {risk_level}\n"
    report += f"  {interpretation}\n"
    
    report += f"""
📋 RECOMMENDATIONS
───────────────────────────────────────────────────────────────────────────────────────
  1. Validate high-risk off-targets experimentally
  2. Consider gRNA design modifications if risk is high
  3. Monitor for off-target effects in tumor samples
  4. Use complementary prediction methods for validation

⚙️  TECHNICAL DETAILS
───────────────────────────────────────────────────────────────────────────────────────
  Genome size assumed:     3.1 × 10⁹ bp
  PAM sequence:           NGG (N = any nucleotide)
  Mismatch tolerance:     Up to 4 mismatches analyzed
  Algorithm:              Exhaustive sequence scanning
"""
    
    return report

def main():
    """Main pipeline function."""
    parser = argparse.ArgumentParser(
        description="Simple CancerPAM: Tumor-specific off-target risk analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic analysis with sample chromosome (chr22)
  python simple_cancerpam.py --gRNA CCGCCAGCGCCGTCTACGTG --tumor_snvs 60000
  
  # Complete genome analysis
  python simple_cancerpam.py --gRNA CCGCCAGCGCCGTCTACGTG --tumor_snvs 60000 --full_genome
  
  # Save results to file
  python simple_cancerpam.py --gRNA CCGCCAGCGCCGTCTACGTG --tumor_snvs 60000 --full_genome --output results.txt
  
  # Use custom genome file
  python simple_cancerpam.py --gRNA CCGCCAGCGCCGTCTACGTG --tumor_snvs 60000 --genome /path/to/genome.fa
        """
    )
    
    parser.add_argument("--gRNA", required=True, 
                       help="20-nucleotide guide RNA sequence (ACGT only)")
    parser.add_argument("--tumor_snvs", type=int, required=True,
                       help="Number of tumor-specific SNVs")
    parser.add_argument("--max_mismatches", type=int, default=4,
                       help="Maximum mismatches to analyze (default: 4)")
    parser.add_argument("--genome", type=str,
                       help="Path to reference genome FASTA file")
    parser.add_argument("--full_genome", action="store_true",
                       help="Use complete human genome (hg38) instead of sample chromosome")
    parser.add_argument("--output", type=str,
                       help="Output file for results")
    parser.add_argument("--no_multiprocessing", action="store_true",
                       help="Disable multiprocessing (useful for debugging)")
    
    args = parser.parse_args()
    
    # Validate inputs
    grna = args.gRNA.upper()
    if not validate_grna(grna):
        sys.exit(1)
    
    if args.tumor_snvs <= 0:
        print("❌ Number of tumor SNVs must be positive")
        sys.exit(1)
    
    # Set up genome file
    if args.genome:
        genome_file = args.genome
        print(f"📁 Using custom genome file: {genome_file}")
    elif args.full_genome:
        # Use complete human genome
        genome_file = Path("CancerPAM_tools") / "hg38.fa.gz"
        genome_file.parent.mkdir(exist_ok=True)
        print("🧬 Using complete human genome (hg38)")
        print("   ⚠️  This will download ~900 MB and analyze ~3.1 GB of sequence data")
        print("   ⏱️  Analysis may take 30-60 minutes depending on your system")
        
        # Confirm with user
        if not genome_file.exists():
            confirm = input("\nProceed with complete genome download and analysis? (y/N): ").strip().lower()
            if confirm != 'y':
                print("❌ Analysis cancelled by user")
                sys.exit(0)
    else:
        # Use default path or download sample
        genome_file = Path("CancerPAM_tools") / "chr22.fa.gz"
        if not genome_file.exists():
            genome_file.parent.mkdir(exist_ok=True)
        print("📊 Using sample chromosome (chr22) for analysis")
        print("   💡 Use --full_genome flag to analyze the complete human genome")
    
    # Initialize analyzer
    try:
        analyzer = SimplePAMAnalyzer(genome_file)
    except Exception as e:
        print(f"❌ Failed to initialize analyzer: {e}")
        sys.exit(1)
    
    # Run analysis
    start_time = time.time()
    try:
        print("\n🚀 Starting off-target analysis...")
        print(f"   Start time: {time.strftime('%Y-%m-%d %H:%M:%S')}")
        
        pam_counts, no_pam_sites = analyzer.analyze_off_targets(
            grna, 
            max_mismatches=args.max_mismatches,
            use_multiprocessing=not args.no_multiprocessing
        )
        
        analysis_time = time.time() - start_time
        print(f"⏱️  Analysis completed in {analysis_time:.1f} seconds")
        print(f"   End time: {time.strftime('%Y-%m-%d %H:%M:%S')}")
        
        # Show detailed results
        print("\n📊 Detailed analysis results:")
        print(f"   PAM-adjacent sites found: {sum(pam_counts.values())}")
        for mismatches in sorted(pam_counts.keys()):
            print(f"     {mismatches} mismatches: {pam_counts[mismatches]} sites")
        
        print(f"   No-PAM sites found: {len(no_pam_sites)}")
        
        # Show some examples of no-PAM sites
        if no_pam_sites:
            print("   Examples of no-PAM sites:")
            for i, (chrom, pos, strand, mm) in enumerate(no_pam_sites[:5]):
                print(f"     {i+1}. {chrom}:{pos} ({strand}) - {mm} mismatches")
        
    except Exception as e:
        print(f"❌ Analysis failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
    
    # Calculate risks
    print("\n🧮 Calculating risks...")
    n_offtargets_4 = pam_counts.get(4, 0)
    no_pam_count = len(no_pam_sites)
    
    print(f"   Input for rescue risk calculation:")
    print(f"     4-mismatch sites: {n_offtargets_4}")
    print(f"     Tumor SNVs: {args.tumor_snvs}")
    
    rescue_lambda, rescue_prob = calculate_rescue_risk(n_offtargets_4, args.tumor_snvs)
    print(f"   Rescue risk: λ = {rescue_lambda:.6f}, P = {rescue_prob:.4%}")
    
    print(f"   Input for novel PAM risk calculation:")
    print(f"     No-PAM sites: {no_pam_count}")
    print(f"     Tumor SNVs: {args.tumor_snvs}")
    
    novel_pam_lambda, novel_pam_prob = calculate_novel_pam_risk(no_pam_count, args.tumor_snvs)
    print(f"   Novel PAM risk: λ = {novel_pam_lambda:.6f}, P = {novel_pam_prob:.4%}")
    
    combined_risk = rescue_prob + novel_pam_prob
    print(f"   Combined risk: {combined_risk:.4%}")
    
    # Generate and display report
    print("\n📝 Generating report...")
    report = generate_detailed_report(
        grna, args.tumor_snvs, pam_counts, no_pam_count,
        rescue_lambda, rescue_prob, novel_pam_lambda, novel_pam_prob
    )
    
    print(report)
    
    # Save to file if requested
    if args.output:
        with open(args.output, 'w') as f:
            f.write(report)
        print(f"📄 Report saved to {args.output}")

if __name__ == "__main__":
    main()
