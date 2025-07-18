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

def build_bowtie_index(bowtie_dir, fasta_path, index_prefix):
    os_platform = platform.system()
    bowtie_build_name = "bowtie-build.exe" if os_platform == "Windows" else "bowtie-build"

    # Find the actual extracted folder (e.g., bowtie-1.2.3-macos-x86_64) inside bowtie_dir
    extracted_bowtie_folder = None
    for item in bowtie_dir.iterdir():
        if item.is_dir() and "bowtie" in item.name and "1.2.3" in item.name:
            extracted_bowtie_folder = item
            break
    
    if not extracted_bowtie_folder:
        print(f"Error: Could not find extracted Bowtie folder inside {bowtie_dir}")
        sys.exit(1)

    bowtie_build_path = find_executable_in_dir(extracted_bowtie_folder, [bowtie_build_name])
    
    if not bowtie_build_path:
        raise FileNotFoundError(f"{bowtie_build_name} not found in {extracted_bowtie_folder} or not executable.")
    
    cmd = [str(bowtie_build_path), str(fasta_path), str(index_prefix)]
    print(f"Building Bowtie index with command: {' '.join(cmd)}")
    try:
        # Run from the extracted folder to ensure relative paths work if needed by bowtie-build
        subprocess.run(cmd, check=True, cwd=str(extracted_bowtie_folder))
        print("Bowtie index built successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error building Bowtie index: {e}")
        sys.exit(1)

def build_bwa_index(bwa_dir, fasta_path):
    os_platform = platform.system()
    
    bwa_exec_name = "bwa.exe" if os_platform == "Windows" else "bwa"
    bwa_exec_path = find_executable_in_dir(bwa_dir, [bwa_exec_name])
    
    if not bwa_exec_path: # If not found, try compiling if not Windows
        if os_platform == "Windows":
            print("\nNOTE: BWA automatic build is not supported on Windows.")
            print("Please install BWA manually (e.g., via WSL/MinGW) and ensure 'bwa.exe' is in your system PATH or place it in the 'CancerPAM_tools/bwa' directory.")
            return None # Indicate failure or skipped
        else:
            print(f"BWA executable not found in {bwa_dir}. Attempting to compile BWA...")
            # Assume bwa source is extracted into a folder like 'bwa-0.7.17'
            bwa_source_folder = None
            for item in bwa_dir.iterdir():
                if item.is_dir() and "bwa" in item.name and "0.7.17" in item.name: # Check for the specific version folder
                    bwa_source_folder = item
                    break
            
            if not bwa_source_folder:
                print(f"Error: Could not find extracted BWA source folder inside {bwa_dir}")
                sys.exit(1)

            try:
                print(f"Navigating to {bwa_source_folder} and running 'make'...")
                subprocess.run(["make"], cwd=str(bwa_source_folder), check=True)
                bwa_exec_path = find_executable_in_dir(bwa_source_folder, [bwa_exec_name]) # Find compiled executable
                if not bwa_exec_path:
                    raise FileNotFoundError(f"BWA executable not found after compilation in {bwa_source_folder}")
                print("BWA compiled successfully.")
            except subprocess.CalledProcessError as e:
                print(f"Error compiling BWA: {e}")
                print("Please try to compile BWA manually from its source code (in 'CancerPAM_tools/bwa').")
                return None # Indicate compilation failure
    
    if bwa_exec_path: # Only proceed if BWA executable is found/compiled
        cmd = [str(bwa_exec_path), "index", str(fasta_path)]
        print(f"Building BWA index with command: {' '.join(cmd)}")
        try:
            subprocess.run(cmd, check=True)
            print("BWA index built successfully.")
            return bwa_exec_path # Return path to the executable
        except subprocess.CalledProcessError as e:
            print(f"Error building BWA index: {e}")
            return None # Indicate index building failure
    return None # BWA executable not available

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

    # Setup options
    print("\n🛠️  Setup Options:")
    print("1. Simple Python-only version (recommended)")
    print("2. Enhanced version with alignment tools")
    print("3. Full setup with all tools")
    
    choice = input("\nChoose setup option (1-3) [1]: ").strip() or "1"
    
    if choice == "1":
        setup_simple_version(base_dir)
    elif choice == "2":
        setup_enhanced_version(base_dir)
    elif choice == "3":
        setup_full_version(base_dir)
    else:
        print("Invalid choice. Using simple version.")
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
    print("python simple_cancerpam.py --gRNA CCGCCAGCGCCGTCTACGTG --tumor_snvs 60000")

def setup_enhanced_version(base_dir):
    """Setup enhanced version with alignment tools."""
    print("\n🚀 Setting up Enhanced version with alignment tools...")
    
    # Setup simple version first
    setup_simple_version(base_dir)
    
    # Try to setup BWA
    print("\n🔧 Setting up BWA...")
    bwa_success = setup_bwa_enhanced(base_dir)
    
    # Try to setup Bowtie2
    print("\n🔧 Setting up Bowtie2...")
    bowtie2_success = setup_bowtie2_enhanced(base_dir)
    
    print(f"\n✅ Enhanced setup complete!")
    print(f"BWA available: {'✓' if bwa_success else '✗'}")
    print(f"Bowtie2 available: {'✓' if bowtie2_success else '✗'}")

def setup_full_version(base_dir):
    """Setup full version with all tools."""
    print("\n🚀 Setting up Full version with all tools...")
    
    # Setup enhanced version first
    setup_enhanced_version(base_dir)
    
    # Cas-OFFinder instructions
    print("\n🔧 Cas-OFFinder Setup (Manual):")
    print("For maximum accuracy, you can optionally install Cas-OFFinder:")
    print("1. Install OpenCL drivers for your system")
    print("2. Download from: http://www.rgenome.net/cas-offinder/")
    print("3. Place executable in CancerPAM_tools/cas-offinder/")
    
    casoffinder_dir = base_dir / "cas-offinder"
    casoffinder_dir.mkdir(exist_ok=True)
    
    print(f"✅ Full setup complete!")
    print(f"📁 Cas-OFFinder directory created: {casoffinder_dir}")

def setup_bwa_enhanced(base_dir):
    """Setup BWA if possible."""
    # Check if BWA is in PATH
    if shutil.which("bwa"):
        print("✓ BWA found in system PATH")
        return True
    
    # Try to download and compile BWA
    bwa_dir = base_dir / "bwa"
    bwa_dir.mkdir(exist_ok=True)
    
    bwa_url = "https://github.com/lh3/bwa/archive/refs/tags/v0.7.17.tar.gz"
    bwa_archive = base_dir / "bwa.tar.gz"
    
    try:
        if not bwa_archive.exists():
            print("⬇️  Downloading BWA...")
            urllib.request.urlretrieve(bwa_url, bwa_archive)
        
        print("📦 Extracting BWA...")
        with tarfile.open(bwa_archive, 'r:gz') as tar:
            tar.extractall(bwa_dir)
        
        # Find source directory
        bwa_src = None
        for item in bwa_dir.iterdir():
            if item.is_dir() and "bwa" in item.name:
                bwa_src = item
                break
        
        if bwa_src and platform.system() != "Windows":
            print("🔨 Compiling BWA...")
            subprocess.run(["make"], cwd=bwa_src, check=True, 
                         stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            
            bwa_exec = bwa_src / "bwa"
            if bwa_exec.exists():
                print("✅ BWA compiled successfully")
                return True
    
    except Exception as e:
        print(f"⚠️  BWA setup failed: {e}")
    
    return False

def setup_bowtie2_enhanced(base_dir):
    """Setup Bowtie2 if possible."""
    # Check if Bowtie2 is in PATH
    if shutil.which("bowtie2"):
        print("✓ Bowtie2 found in system PATH")
        return True
    
    # Try to download precompiled Bowtie2
    bowtie2_dir = base_dir / "bowtie2"
    bowtie2_dir.mkdir(exist_ok=True)
    
    os_platform = platform.system()
    if os_platform == "Darwin":
        bowtie2_url = "https://github.com/BenLangmead/bowtie2/releases/download/v2.5.1/bowtie2-2.5.1-macos-x86_64.zip"
    elif os_platform == "Linux":
        bowtie2_url = "https://github.com/BenLangmead/bowtie2/releases/download/v2.5.1/bowtie2-2.5.1-linux-x86_64.zip"
    else:
        print("⚠️  Bowtie2 auto-setup not supported on Windows")
        return False
    
    try:
        bowtie2_archive = base_dir / "bowtie2.zip"
        if not bowtie2_archive.exists():
            print("⬇️  Downloading Bowtie2...")
            urllib.request.urlretrieve(bowtie2_url, bowtie2_archive)
        
        print("📦 Extracting Bowtie2...")
        with zipfile.ZipFile(bowtie2_archive, 'r') as zip_ref:
            zip_ref.extractall(bowtie2_dir)
        
        # Find executable
        for item in bowtie2_dir.rglob("bowtie2"):
            if item.is_file():
                print("✅ Bowtie2 installed successfully")
                return True
    
    except Exception as e:
        print(f"⚠️  Bowtie2 setup failed: {e}")
    
    return False

if __name__ == "__main__":
    main_setup()