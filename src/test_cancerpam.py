#!/usr/bin/env python3
"""
Test script for CancerPAM pipeline
==================================

This script tests the various components of the CancerPAM pipeline
to ensure everything is working correctly.
"""

import sys
import os
from pathlib import Path

def test_imports():
    """Test if all required packages are available."""
    print("🧪 Testing package imports...")
    
    try:
        import numpy
        print("✓ NumPy available")
    except ImportError:
        print("❌ NumPy not available")
        return False
    
    try:
        import tqdm
        print("✓ tqdm available")
    except ImportError:
        print("❌ tqdm not available")
        return False
    
    try:
        import pyfaidx
        print("✓ pyfaidx available")
    except ImportError:
        print("❌ pyfaidx not available")
        return False
    
    try:
        from Bio import SeqIO
        print("✓ BioPython available")
    except ImportError:
        print("❌ BioPython not available")
        return False
    
    return True

def test_simple_pipeline():
    """Test the simple pipeline with a mock sequence."""
    print("\n🧪 Testing simple pipeline...")
    
    # Create a minimal test sequence
    test_dir = Path("test_data")
    test_dir.mkdir(exist_ok=True)
    
    test_fasta = test_dir / "test.fa"
    
    # Create a simple test sequence with known off-targets
    test_sequence = (
        ">chr_test\n"
        "ATCGATCGATCGATCGATCGCCGCCAGCGCCGTCTACGTGGGATCGATCGATCGATCG\n"
        "ATCGATCGATCGATCGATCGCCGCCAGCGCCGTCTACGTCGGATCGATCGATCGATCG\n"
        "ATCGATCGATCGATCGATCGCCGCCAGCGCCGTCTACGTGGGAATCGATCGATCGATCG\n"
    )
    
    with open(test_fasta, 'w') as f:
        f.write(test_sequence)
    
    # Try to import and run the simple analyzer
    try:
        sys.path.insert(0, str(Path.cwd()))
        from simple_cancerpam import SimplePAMAnalyzer
        
        analyzer = SimplePAMAnalyzer(test_fasta)
        
        # Test with a gRNA that should have matches
        test_grna = "CCGCCAGCGCCGTCTACGTG"
        pam_counts, no_pam_sites = analyzer.analyze_off_targets(test_grna, max_mismatches=2)
        
        print(f"✓ Analysis completed")
        print(f"  PAM-adjacent sites: {sum(pam_counts.values())}")
        print(f"  No-PAM sites: {len(no_pam_sites)}")
        
        return True
        
    except Exception as e:
        print(f"❌ Pipeline test failed: {e}")
        return False

def test_risk_calculations():
    """Test risk calculation functions."""
    print("\n🧪 Testing risk calculations...")
    
    try:
        sys.path.insert(0, str(Path.cwd()))
        from simple_cancerpam import calculate_rescue_risk, calculate_novel_pam_risk
        
        # Test rescue risk
        rescue_lambda, rescue_prob = calculate_rescue_risk(
            n_offtargets_4=100, 
            n_mutations=60000
        )
        print(f"✓ Rescue risk calculation: λ={rescue_lambda:.6f}, P={rescue_prob:.4f}")
        
        # Test novel PAM risk
        pam_lambda, pam_prob = calculate_novel_pam_risk(
            no_pam_count=1000,
            n_mutations=60000
        )
        print(f"✓ Novel PAM risk calculation: λ={pam_lambda:.6f}, P={pam_prob:.4f}")
        
        return True
        
    except Exception as e:
        print(f"❌ Risk calculation test failed: {e}")
        return False

def test_sequence_functions():
    """Test basic sequence manipulation functions."""
    print("\n🧪 Testing sequence functions...")
    
    try:
        sys.path.insert(0, str(Path.cwd()))
        from simple_cancerpam import SimplePAMAnalyzer
        
        # Create a dummy analyzer to test functions
        test_fasta = Path("test_data/test.fa")
        if not test_fasta.exists():
            return False
        
        analyzer = SimplePAMAnalyzer(test_fasta)
        
        # Test reverse complement
        seq = "ATCG"
        rc = analyzer.reverse_complement(seq)
        expected = "CGAT"
        
        if rc == expected:
            print(f"✓ Reverse complement: {seq} -> {rc}")
        else:
            print(f"❌ Reverse complement failed: expected {expected}, got {rc}")
            return False
        
        # Test mismatch counting
        mismatches = analyzer.count_mismatches("ATCG", "ATGG")
        if mismatches == 1:
            print(f"✓ Mismatch counting: 1 mismatch detected")
        else:
            print(f"❌ Mismatch counting failed: expected 1, got {mismatches}")
            return False
        
        return True
        
    except Exception as e:
        print(f"❌ Sequence function test failed: {e}")
        return False

def main():
    """Run all tests."""
    print("🧬 CancerPAM Pipeline Test Suite")
    print("=" * 40)
    
    tests = [
        test_imports,
        test_sequence_functions,
        test_risk_calculations,
        test_simple_pipeline
    ]
    
    passed = 0
    failed = 0
    
    for test in tests:
        try:
            if test():
                passed += 1
            else:
                failed += 1
        except Exception as e:
            print(f"❌ Test {test.__name__} crashed: {e}")
            failed += 1
    
    print(f"\n📊 Test Results:")
    print(f"  Passed: {passed}")
    print(f"  Failed: {failed}")
    
    if failed == 0:
        print("🎉 All tests passed! The pipeline is ready to use.")
        print("\nTo run the pipeline:")
        print("python simple_cancerpam.py --gRNA CCGCCAGCGCCGTCTACGTG --tumor_snvs 60000")
    else:
        print("❌ Some tests failed. Please check the setup.")
        return 1
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
