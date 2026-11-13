"""Fail-closed numerical artifacts for the working RDT repository.
Single writer only. Interrupted writes are rejected rather than silently repaired.
Standalone releases must ship this module with their entry points.
"""
from pathlib import Path
import csv
import hashlib
import json
import math
import platform
import sys
import numpy as np
import scipy


def positive(**values):
    if any(not math.isfinite(v) or v <= 0 for v in values.values()):
        raise ValueError(f'Expected finite positive values: {values}')


def check_problem(p, rho_c, rho0, xi0, xi_max, rtol, atol, max_step):
    positive(rho_c=rho_c, rho0=rho0, xi0=xi0, xi_max=xi_max,
             rtol=rtol, atol=atol, max_step=max_step, alpha=p.alpha)
    if not math.isfinite(p.A) or not 0 <= p.A < 1 or xi0 >= xi_max:
        raise ValueError('Require 0 <= A < 1 and 0 < xi0 < xi_max')


def check_surface(sol):
    if not sol.success:
        raise RuntimeError(f'ODE solver failed: {sol.message}')
    if sol.status != 1 or len(sol.t_events[0]) != 1:
        raise RuntimeError('No unique terminal stellar surface found')
    if (not np.isfinite(sol.t).all() or not np.isfinite(sol.y).all()
            or not np.isfinite(sol.y_events[0]).all()):
        raise RuntimeError('Nonfinite stellar solution')
    if (len(sol.t) < 2 or np.any(np.diff(sol.t) <= 0)
            or abs(sol.t[-1]-sol.t_events[0][0]) > 1e-10
            or abs(sol.y[0,-1]) > 1e-8 or sol.y[1,-1] >= 0):
        raise RuntimeError('Invalid surface endpoint')


def run_gates(gates):
    results=[]
    for gate in gates:
        try:
            ok=bool(gate())
        except Exception as exc:
            print(f'{gate.__name__}: FAIL ({type(exc).__name__}: {exc})')
            ok=False
        results.append(ok)
    if not all(results):
        raise RuntimeError('Validation failed; numerical output refused')
    print('ALL GATES PASS')
    return True


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def fingerprint(directory, files, settings):
    directory=Path(directory)
    return {'schema':1, 'source_sha256':{f:sha(directory/f) for f in files},
            'guardrails_sha256':sha(__file__), 'settings':settings,
            'runtime':{'python':sys.version,'numpy':np.__version__,
                       'scipy':scipy.__version__,'platform':platform.platform()}}


def finite_row(row, fields):
    if set(row) != set(fields):
        raise ValueError('Unexpected row columns')
    out={k:float(row[k]) for k in fields}
    if not all(math.isfinite(v) for v in out.values()):
        raise ValueError('Nonfinite output row')
    return out


class Artifact:
    def __init__(self, path, contract, fields, keys, expected, validate):
        self.path=Path(path)
        self.meta=self.path.with_suffix(self.path.suffix+'.manifest.json')
        self.log=self.path.with_suffix(self.path.suffix+'.failures.jsonl')
        self.contract=contract; self.fields=fields; self.keys=keys
        self.expected={tuple(map(float,k)) for k in expected}; self.validate=validate
        if not self.expected or len(self.expected)!=len(expected):
            raise ValueError('Empty or duplicate expected grid')
        self.rows=[]; self.state=None

    def load(self, complete=False):
        if not self.path.exists() and not self.meta.exists():
            if complete: raise ValueError('Missing numerical artifact')
            return []
        if not self.path.exists() or not self.meta.exists():
            raise ValueError('Unmanifested/incomplete artifact; use a new output path')
        data=json.loads(self.meta.read_text())
        if data['contract'] != self.contract:
            raise ValueError('Incompatible source/settings/runtime provenance')
        if data['status'] not in ('partial','failed','complete'):
            raise ValueError('Invalid artifact status')
        if data['sha256'] != sha(self.path):
            raise ValueError('Artifact checksum mismatch (possibly interrupted write)')
        with self.path.open() as f:
            reader=csv.DictReader(f)
            if reader.fieldnames != self.fields: raise ValueError('Invalid CSV header')
            self.rows=[self.checked(r) for r in reader]
        keys=[tuple(r[k] for k in self.keys) for r in self.rows]
        if len(keys)!=len(set(keys)) or len(keys)!=data['rows']:
            raise ValueError('Duplicate rows or count mismatch')
        if data['status']=='complete' and set(keys)!=self.expected:
            raise ValueError('Incomplete grid labeled complete')
        if complete and data['status']!='complete':
            raise ValueError('Partial or failed run cannot be analyzed')
        self.state=data['status']
        return self.rows

    def checked(self, row):
        row=finite_row(row,self.fields)
        if tuple(row[k] for k in self.keys) not in self.expected:
            raise ValueError('Row outside declared parameter grid')
        self.validate(row)
        return row

    def save(self, status):
        if status=='complete' and {tuple(r[k] for k in self.keys) for r in self.rows}!=self.expected:
            raise ValueError('Cannot complete an incomplete grid')
        self.path.parent.mkdir(parents=True,exist_ok=True)
        temp=self.path.with_suffix(self.path.suffix+'.tmp')
        with temp.open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=self.fields);w.writeheader();w.writerows(self.rows)
        temp.replace(self.path)
        metadata={'contract':self.contract,'status':status,'rows':len(self.rows),'sha256':sha(self.path)}
        tmp=self.meta.with_suffix('.tmp');tmp.write_text(json.dumps(metadata,indent=2)+'\n');tmp.replace(self.meta)
        self.state=status

    def append(self, row):
        row=self.checked(row)
        key=tuple(row[k] for k in self.keys)
        if key in {tuple(r[k] for k in self.keys) for r in self.rows}:
            raise ValueError('Duplicate new row')
        self.rows.append(row);self.save('partial')

    def fail(self, point, exc):
        self.log.parent.mkdir(parents=True,exist_ok=True)
        with self.log.open('a') as f:
            f.write(json.dumps({'point':point,'error':str(exc),'type':type(exc).__name__})+'\n')
        self.save('failed')
