"""Supported entry point: require targeted numerical/inference gates before a run."""
from pathlib import Path
import subprocess,sys
ROOT=Path(__file__).resolve().parents[2]
if __name__=='__main__':
    subprocess.run([sys.executable,str(ROOT/'tests/test_step11_inference.py')],check=True)
    subprocess.run([sys.executable,str(Path(__file__).with_name('run.py')),*sys.argv[1:]],check=True)
