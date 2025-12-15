import os
import subprocess
import platform

def run_command(command):
    """Execute a shell command and return the result."""
    try:
        result = subprocess.run(
            command, 
            shell=True, 
            check=True, 
            stdout=subprocess.PIPE, 
            stderr=subprocess.PIPE,
            text=True
        )
        return True, result.stdout
    except subprocess.CalledProcessError as e:
        return False, e.stderr

def is_installed(package_name, check_command):
    """Check if a package is installed."""
    success, _ = run_command(check_command)
    return success

def install_packages():
    """Install make, GSL, and GNU C++ libraries if not present."""
    
    # Check OS
    if platform.system() != "Linux":
        print("This script is designed for Linux systems only.")
        return
    
    # Detect package manager
    if is_installed("apt", "which apt"):
        pkg_manager = "apt"
        install_cmd = "sudo apt update && sudo apt install -y"
        packages = {
            "make": ("which make", "make"),
            "gsl": ("ldconfig -p | grep libgsl", "libgsl-dev"),
            "g++": ("which g++", "g++")
        }
    elif is_installed("yum", "which yum"):
        pkg_manager = "yum"
        install_cmd = "sudo yum install -y"
        packages = {
            "make": ("which make", "make"),
            "gsl": ("ldconfig -p | grep libgsl", "gsl-devel"),
            "g++": ("which g++", "gcc-c++")
        }
    elif is_installed("dnf", "which dnf"):
        pkg_manager = "dnf"
        install_cmd = "sudo dnf install -y"
        packages = {
            "make": ("which make", "make"),
            "gsl": ("ldconfig -p | grep libgsl", "gsl-devel"),
            "g++": ("which g++", "gcc-c++")
        }
    else:
        print("Unsupported package manager. Please install manually.")
        return
    
    print(f"Detected package manager: {pkg_manager}\n")
    
    # Check and install each package
    for name, (check_cmd, package) in packages.items():
        print(f"Checking {name}...")
        if is_installed(name, check_cmd):
            print(f"✓ {name} is already installed.\n")
        else:
            print(f"✗ {name} not found. Installing {package}...")
            success, output = run_command(f"{install_cmd} {package}")
            if success:
                print(f"✓ {package} installed successfully.\n")
            else:
                print(f"✗ Failed to install {package}. Error: {output}\n")

if __name__ == "__main__":
    print("=== C++ Development Tools Installer ===\n")
    install_packages()
    print("Installation check complete!")

    current_dir = os.getcwd()
    print(current_dir)
    os.chdir("06receptors/base/")
    print(os.getcwd())
    
    os.system("make all")
    os.system("make run")
    
    os.chdir(current_dir+'/24receptors/base/')
    print(os.getcwd())
    
    os.system("make all")
    os.system("make run")
    
    
    