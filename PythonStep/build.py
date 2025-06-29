#!/usr/bin/env python3
"""
Build script for PythonStep package
"""

import os
import sys
import subprocess
import shutil
from pathlib import Path


def check_dependencies():
    """Check if required dependencies are available"""
    print("Checking dependencies...")
    
    # Check Python version
    if sys.version_info < (3, 7):
        print("❌ Python 3.7+ required")
        return False
    print(f"✅ Python {sys.version_info.major}.{sys.version_info.minor}")
    
    # Check required Python packages
    required_packages = ['numpy', 'cython', 'setuptools']
    for package in required_packages:
        try:
            __import__(package)
            print(f"✅ {package}")
        except ImportError:
            print(f"❌ {package} not found")
            return False
    
    # Check C++ compiler
    try:
        result = subprocess.run(['g++', '--version'], 
                              capture_output=True, text=True)
        if result.returncode == 0:
            print("✅ g++ compiler found")
        else:
            print("❌ g++ compiler not found")
            return False
    except FileNotFoundError:
        print("❌ g++ compiler not found")
        return False
    
    return True


def check_computree_libraries():
    """Check if COMPUTREE libraries are available"""
    print("\nChecking COMPUTREE libraries...")
    
    # Common library paths
    lib_paths = [
        '/usr/local/lib',
        '/usr/lib',
        '/opt/computree/lib',
        '/usr/local/computree/lib',
    ]
    
    # Required libraries
    required_libs = [
        'libctlibplugin.so',
        'libctlibclouds.so',
        'libctlibstep.so',
        'libctlibstepaddon.so',
        'libctlibio.so',
        'libctlibfilters.so',
        'libctlibaction.so',
        'libctlibstdactions.so',
        'libctlibmath.so',
        'libctliblas.so',
    ]
    
    found_libs = []
    for lib in required_libs:
        found = False
        for lib_path in lib_paths:
            lib_file = os.path.join(lib_path, lib)
            if os.path.exists(lib_file):
                found = True
                found_libs.append(lib)
                print(f"✅ {lib} found in {lib_path}")
                break
        if not found:
            print(f"❌ {lib} not found")
    
    if len(found_libs) == len(required_libs):
        print("✅ All COMPUTREE libraries found")
        return True
    else:
        print(f"⚠️  Only {len(found_libs)}/{len(required_libs)} libraries found")
        print("You may need to update library paths in setup.py")
        return False


def update_setup_paths():
    """Update setup.py with detected library paths"""
    print("\nUpdating setup.py with detected paths...")
    
    setup_file = "setup.py"
    if not os.path.exists(setup_file):
        print("❌ setup.py not found")
        return False
    
    # Read current setup.py
    with open(setup_file, 'r') as f:
        content = f.read()
    
    # Find library paths (this is a simple approach)
    # In a real implementation, you might want to use a more robust method
    print("⚠️  Please manually update library paths in setup.py if needed")
    print("Common paths to check:")
    print("  - /usr/local/lib")
    print("  - /usr/lib")
    print("  - /opt/computree/lib")
    print("  - /usr/local/computree/lib")
    
    return True


def build_package():
    """Build the PythonStep package"""
    print("\nBuilding PythonStep package...")
    
    try:
        # Clean previous builds
        if os.path.exists("build"):
            shutil.rmtree("build")
        if os.path.exists("PythonStep.egg-info"):
            shutil.rmtree("PythonStep.egg-info")
        
        # Build extension
        result = subprocess.run([
            sys.executable, "setup.py", "build_ext", "--inplace"
        ], capture_output=True, text=True)
        
        if result.returncode == 0:
            print("✅ Package built successfully")
            return True
        else:
            print("❌ Build failed")
            print("Error output:")
            print(result.stderr)
            return False
            
    except Exception as e:
        print(f"❌ Build error: {e}")
        return False


def install_package():
    """Install the PythonStep package"""
    print("\nInstalling PythonStep package...")
    
    try:
        result = subprocess.run([
            sys.executable, "-m", "pip", "install", "-e", "."
        ], capture_output=True, text=True)
        
        if result.returncode == 0:
            print("✅ Package installed successfully")
            return True
        else:
            print("❌ Installation failed")
            print("Error output:")
            print(result.stderr)
            return False
            
    except Exception as e:
        print(f"❌ Installation error: {e}")
        return False


def run_tests():
    """Run the test suite"""
    print("\nRunning tests...")
    
    try:
        result = subprocess.run([
            sys.executable, "-m", "unittest", "discover", "tests"
        ], capture_output=True, text=True)
        
        if result.returncode == 0:
            print("✅ Tests passed")
            return True
        else:
            print("❌ Tests failed")
            print("Test output:")
            print(result.stdout)
            print("Error output:")
            print(result.stderr)
            return False
            
    except Exception as e:
        print(f"❌ Test error: {e}")
        return False


def main():
    """Main build function"""
    print("=== PythonStep Build Script ===\n")
    
    # Check dependencies
    if not check_dependencies():
        print("\n❌ Dependency check failed. Please install missing dependencies.")
        return 1
    
    # Check COMPUTREE libraries
    check_computree_libraries()
    
    # Update setup paths
    update_setup_paths()
    
    # Build package
    if not build_package():
        print("\n❌ Build failed. Please check the error messages above.")
        return 1
    
    # Install package
    if not install_package():
        print("\n❌ Installation failed. Please check the error messages above.")
        return 1
    
    # Run tests
    if not run_tests():
        print("\n❌ Tests failed. Please check the error messages above.")
        return 1
    
    print("\n🎉 PythonStep package built and installed successfully!")
    print("\nYou can now use it in Python:")
    print("  import PythonStep as ps")
    print("  grid = ps.create_grid_from_txt('twoCylinders.txt')")
    
    return 0


if __name__ == "__main__":
    sys.exit(main()) 