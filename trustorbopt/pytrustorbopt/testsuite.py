import subprocess
import sys
from importlib import resources

with resources.path("pytrustorbopt", "testsuite") as exe_path:
    if not exe_path.exists():
        print("Error: Testsuite executable not found.")
        sys.exit(1)
    subprocess.check_call([str(exe_path)])
