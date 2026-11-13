"""Versioned solar scan. Single writer; explicit new output path required.
Historical unmanifested CSVs are read-only evidence, never resumable checkpoints.
This producer retains legacy measurement scoring; step 11 changes inference.
"""
import argparse
import json
import math
from pathlib import Path
import numpy as np
from model import (Params, calibrate, flux_ratios, chi2_fluxes,
                   sound_speed_deviation, BETA, EXP_R, EXP_SIGMA, THETA_FLOOR)
from numerical_guardrails import Artifact, fingerprint, positive
import validate

HERE=Path(__file__).resolve().parent
A_VALS=np.round(np.linspace(0,.040,41),6)
ALPHA_VALS=np.round(np.geomspace(.02,.50,33),6)
LE_OPTS=dict(rtol=1e-9,atol=1e-11,max_step=.05)
FIELDS=['A','alpha','xi1','rho_c','Tc_ratio','R_pp','R_pep','R_be7','R_b8','chi2','max_abs_dcs']


def check_row(r):
    positive(**{k:r[k] for k in ['xi1','rho_c','Tc_ratio',*[f'R_{b}' for b in BETA]]})
    if r['max_abs_dcs']<0: raise ValueError('Negative sound-speed maximum')
    score=chi2_fluxes({b:r[f'R_{b}'] for b in BETA})
    if not math.isclose(score,r['chi2'],rel_tol=1e-12,abs_tol=1e-12):
        raise ValueError('Cached score does not match its source ratios and current inputs')
    if r['A']==0 and (any(r[f'R_{b}']!=1 for b in BETA) or r['Tc_ratio']!=1 or r['max_abs_dcs']!=0):
        raise ValueError('Invalid null row')


def artifact(path, a_values, alpha_values):
    aa=list(map(float,a_values));al=list(map(float,alpha_values))
    if (not aa or not al or aa!=sorted(set(aa)) or al!=sorted(set(al))
            or any(not math.isfinite(a) or not 0<=a<=.04 for a in aa)
            or any(not math.isfinite(a) or not .02<=a<=.5 for a in al) or aa[0]!=0):
        raise ValueError('Invalid canonical grid (requires A=0)')
    settings={'A':aa,'alpha':al,'solver':LE_OPTS,'theta_floor':THETA_FLOOR,
              'reference':'a','calibration_tol':1e-12,'max_iter':60,
              'score':'legacy_measurement_v1','observed':EXP_R,'sigma':EXP_SIGMA,
              'beta':BETA,'sound_speed_grid':[.01,.98,400]}
    contract=fingerprint(HERE,['model.py','scan.py','validate.py'],settings)
    return Artifact(path,contract,FIELDS,['A','alpha'],[(a,l) for a in aa for l in al],check_row)


def load_complete(path):
    meta=Path(path).with_suffix(Path(path).suffix+'.manifest.json')
    if not meta.exists(): raise ValueError('Unmanifested scan cannot enter current analysis')
    settings=json.loads(meta.read_text())['contract']['settings']
    data=artifact(path,settings['A'],settings['alpha'])
    return data.load(complete=True)


def main(output, a_values=A_VALS, alpha_values=ALPHA_VALS):
    data=artifact(output,a_values,alpha_values)
    data.load()  # refuse old/incompatible/corrupt output before any mutation
    point={'phase':'validation'}
    try:
        if validate.run_all() is not True: raise RuntimeError('Validation did not pass')
        point={'phase':'standard_calibration'}
        std=calibrate(Params(0,.1),**LE_OPTS)
        done={(r['A'],r['alpha']) for r in data.rows}
        for A in a_values:
            for alpha in alpha_values:
                if (A,alpha) in done: continue
                point={'A':float(A),'alpha':float(alpha)}
                if A==0:
                    ratios={b:1. for b in BETA}
                    row=dict(A=A,alpha=alpha,xi1=std.xi1,rho_c=std.rho_c,Tc_ratio=1.,max_abs_dcs=0.)
                else:
                    mod=calibrate(Params(A,alpha),**LE_OPTS)
                    ratios=flux_ratios(mod,std)
                    _,dcs=sound_speed_deviation(mod,std)
                    row=dict(A=A,alpha=alpha,xi1=mod.xi1,rho_c=mod.rho_c,
                             Tc_ratio=mod.Tc_unit/std.Tc_unit,max_abs_dcs=float(np.max(np.abs(dcs))))
                row.update({f'R_{b}':ratios[b] for b in BETA})
                row['chi2']=chi2_fluxes(ratios)
                data.append(row)
        data.save('complete')
    except Exception as exc:
        data.fail(point,exc)
        raise
    print(f'Complete: {len(data.rows)} points; {output}')
    return data.rows


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',required=True,help='New CSV path, or compatible manifested checkpoint')
    args=parser.parse_args()
    main(args.output)
