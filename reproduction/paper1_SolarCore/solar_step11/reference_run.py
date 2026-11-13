"""Explicit tighter-ODE reference-prescription run after recorded failure.
The unchanged runner fingerprints the resulting Settings in its manifest.
"""
from functools import partial
from pathlib import Path
import subprocess,sys
import run
from solver import Settings
if __name__=='__main__':
    subprocess.run([sys.executable,str(Path(__file__).resolve().parents[2]/'tests/test_step11_inference.py')],check=True)
    run.Settings=partial(Settings,rtol=2e-13,atol=2e-15,max_step=.02)
    run.run(sys.argv[1],reference='Rsun',na=41,nalpha=33)
