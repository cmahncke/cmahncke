import subprocess
import sys


# This file checks and installs all required packages from pip.
# This function checks and installs packages.
def install(package):
    subprocess.check_call([sys.executable, "-m", "pip", "install", package])
    return 0


install("numpy")
install("requests")
install("xlsxwriter")

print("\n\nAll required packages installed\nmain.py is ready for use!")