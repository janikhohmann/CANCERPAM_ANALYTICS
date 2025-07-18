import sys
import os
import subprocess
import shutil
import urllib.request
import zipfile
import tarfile
import gzip
from pathlib import Path
import importlib
import platform

def check_python_version():
    if sys.version_info < (3, 7):
        print("Python 3.7 or higher is required. Aborting.")
        sys.exit(1)

def install_python_packages(packages):
    """Install Python packages with better error handling."""
    missing_packages = []
    
    for pkg in packages:
        pkg_name = pkg.split('>=')[0].split('==')[0]  # Handle version specs
        try:
            importlib.import_module(pkg_name)
            print(f"✓ Package '{pkg_name}' already installed.")
        except ImportError:
            missing_packages.append(pkg)
    
    if missing_packages:
        print(f"📦 Installing packages: {', '.join(missing_packages)}")
        try:
            subprocess.check_call([
                sys.executable, "-m", "pip", "install", "--upgrade"
            ] + missing_packages)
            print("✅ All packages installed successfully!")
        except subprocess.CalledProcessError as e:
            print(f"❌ Error installing packages: {e}")
            print("Please try installing manually:")
            print(f"pip install {' '.join(missing_packages)}")
            sys.exit(1)

def download_file(url, dest_path):
    if dest_path.exists():
        print(f"{dest_path} already exists, skipping download.")
        return
    print(f"Downloading {url} ...")
    try:
        urllib.request.urlretrieve(url, dest_path)
        print("Download completed.")
    except urllib.error.URLError as e:
        print(f"Error downloading {url}: {e.reason}. Please check your internet connection or the URL.")
        # Do not exit here for critical components, just inform
        return False # Indicate failure
    return True # Indicate success

def unzip_file(zip_path, extract_dir):
    print(f"Extracting {zip_path} ...")
    try:
        with zipfile.ZipFile(zip_path, 'r') as zip_ref:
            zip_ref.extractall(extract_dir)
        print(f"Extracted to {extract_dir}")
    except zipfile.BadZipFile:
        print(f"Error: {zip_path} is not a valid zip file. It might be corrupted or incomplete.")
        sys.exit(1)

def untar_file(tar_path, extract_dir):
    print(f"Extracting {tar_path} ...")
    try:
        with tarfile.open(tar_path, 'r:bz2') as tar_ref: # For .tar.bz2
            tar_ref.extractall(extract_dir)
        print(f"Extracted to {extract_dir}")
    except tarfile.ReadError:
        print(f"Error: {tar_path} is not a valid tar.bz2 file. It might be corrupted or incomplete.")
        sys.exit(1)

def find_executable_in_dir(directory, possible_names):
    """Searches for an executable in a directory given possible names (e.g., 'bwa', 'bwa.exe')."""
    for name in possible_names:
        exec_path = directory / name
        if exec_path.exists():
            # On Linux/macOS, check if it's executable
            if platform.system() != "Windows" and not os.access(exec_path, os.X_OK):
                try:
                    subprocess.run(["chmod", "+x", str(exec_path)], check=True)
                    print(f"Permissions set for {exec_path}")
                    return exec_path
                except subprocess.CalledProcessError:
                    print(f"Warning: Could not set executable permissions for {exec_path}")
                    continue # Try next name or fail
            return exec_path
    return None


def main_setup():
    """Enhanced setup with multiple options."""
    print("🧬 CancerPAM Tools Setup")
    print("=" * 50)
    
    check_python_version()

    base_dir = Path("CancerPAM_tools")
    base_dir.mkdir(exist_ok=True)

    print("\n📦 Installing Python packages...")
    # Essential packages for the simple version
    essential_packages = [
        'numpy>=1.19.0',
        'tqdm>=4.50.0',
        'pyfaidx>=0.6.0',  # For FASTA handling
        'biopython>=1.78'   # For sequence operations
    ]
    
    # Optional packages for enhanced functionality
    optional_packages = [
        'pysam>=0.16.0',   # For SAM/BAM file handling
        'pandas>=1.1.0',   # For data analysis
        'scipy>=1.5.0'     # For statistical functions
    ]
    
    install_python_packages(essential_packages)
    
    # Try to install optional packages
    print("\n📦 Installing optional packages...")
    for pkg in optional_packages:
        try:
            subprocess.check_call([sys.executable, "-m", "pip", "install", pkg], 
                                 stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            print(f"✓ Installed {pkg}")
        except subprocess.CalledProcessError:
            print(f"⚠️  Optional package {pkg} failed to install (continuing anyway)")

    setup_simple_version(base_dir)

def setup_simple_version(base_dir):
    """Setup simple Python-only version."""
    print("\n🚀 Setting up Simple Python-only version...")
    
    # Download sample genome for testing
    sample_genome = base_dir / "chr22.fa.gz"
    if not sample_genome.exists():
        print("⬇️  Downloading sample genome (chromosome 22)...")
        chr22_url = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr22.fa.gz"
        try:
            urllib.request.urlretrieve(chr22_url, sample_genome)
            print("✅ Sample genome downloaded")
        except Exception as e:
            print(f"⚠️  Failed to download sample genome: {e}")
            print("You can download it manually or use your own genome file")
    
    # Create simple config
    config_content = f'''"""
Simple CancerPAM Configuration
"""

# Sample genome (chromosome 22 for testing)
SAMPLE_GENOME = r"{sample_genome}"

# For full analysis, download complete genome:
# wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
# FULL_GENOME = r"{base_dir}/hg38.fa.gz"

# Method to use
ANALYSIS_METHOD = "python_native"

# Configuration version
CONFIG_VERSION = "1.0_simple"
'''
    
    config_file = base_dir / "simple_config.py"
    with open(config_file, 'w') as f:
        f.write(config_content)
    
    print(f"✅ Simple setup complete!")
    print(f"📁 Tools directory: {base_dir}")
    print(f"📄 Configuration: {config_file}")
    print(f"🧬 Sample genome: {sample_genome}")
    print("\n🚀 You can now run:")
    print("python cancerpam_analytics.py --gRNA CCGCCAGCGCCGTCTACGTG --tumor_snvs 60000")


if __name__ == "__main__":
    main_setup()