"""Paper II audit: independent flux ODE, no imported production physics.
Run: python3 paper2_WhiteDwarfs/audit_2026_11_11/audit.py
Writes only alongside this script. Failures raise/nonzero; no sign gate.
"""
from pathlib import Path
import contextlib, csv, hashlib, io, json, platform, sys
import numpy as np
import scipy
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PROD = HERE.parent / 'wd_v5'
sys.path.insert(0, str(PROD))
import model as production
import validate


def flux(A, alpha, rho, *, rtol=2e-11, atol=2e-13, step=.01, eps=1e-6):
    # Independently transcribed stated n=3/2, mu_e=2, cgs model.
    n, G, K, rho_ref = 1.5, 6.67430e-8, 1.0036e13 / 2**(5/3), 150.
    def delta(t):
        x = (rho * max(t, 0.)**n / rho_ref)**alpha
        return -3*A*x/(1+x)
    dc = delta(1.)
    def rhs(xi, y):
        th = max(y[0], 0.)
        w = xi**(2+delta(th))
        return [-y[1]/w, w*th**n, xi**2*th**n]
    def surface(xi,y): return y[0]
    surface.terminal, surface.direction = True, -1
    initial = [1-eps**2/(2*(3+dc)), eps**(3+dc)/(3+dc), eps**3/3]
    sol = solve_ivp(rhs, (eps, 10.), initial, method='RK45', events=surface,
                    rtol=rtol, atol=atol, max_step=step)
    if not sol.success or len(sol.t_events[0]) != 1 or not np.isfinite(sol.y).all():
        raise RuntimeError(sol.message)
    xi = sol.t_events[0][0]
    th,q,I = sol.y_events[0][0]
    a = np.sqrt((n+1)*K*rho**(1/n-1)/(4*np.pi*G))
    return dict(xi=float(xi), q=float(q), I=float(I),
                M=float(4*np.pi*a**3*rho*I), R=float(a*xi))


def shifts(s,m):
    dr, dm = np.log(m['R']/s['R']), np.log(m['M']/s['M'])
    return np.array([dr, dm, dr+dm/3])


def prod_pair(A, alpha, rho, **kw):
    s,m = production.mr_pair(production.Params(A,alpha),rho,**kw)
    return shifts(vars(s),vars(m))


def run():
    checks = {}
    gate_text = io.StringIO()
    with contextlib.redirect_stdout(gate_text):
        gates = [getattr(validate,f'gate{i}')() for i in range(1,7)]
    (HERE/'production_gates.txt').write_text(gate_text.getvalue())
    checks['production_six_gates'] = bool(all(gates))
    b = flux(0.,.15,1e6)
    checks['baseline_xi'] = abs(b['xi']/3.65375374-1) < 1e-8
    checks['baseline_integral'] = abs(b['I']/2.71406-1) < 3e-6
    checks['baseline_flux_equals_matter_integral'] = abs(b['I']-b['q']) < 1e-11
    baseline = {str(r):flux(0.,.15,r) for r in (1e5,1e6,3e6)}
    slopes = dict(R=float(np.log(baseline['3000000.0']['R']/baseline['100000.0']['R'])/np.log(30)),
                  M=float(np.log(baseline['3000000.0']['M']/baseline['100000.0']['M'])/np.log(30)))
    checks['baseline_scalings'] = abs(slopes['R']+1/6)<1e-10 and abs(slopes['M']-.5)<1e-10

    # All archived grid rows: independent flux solve and production FAST replay.
    rows = list(csv.DictReader((PROD/'results_predictions.csv').open()))
    grid=[]
    for row in rows:
        a,A,rho = (float(row[k]) for k in ('alpha','A','rho_c'))
        std,mod = flux(0,a,rho),flux(A,a,rho)
        independent = shifts(std,mod)
        saved = np.array([float(row[k]) for k in ('dlnR','dlnM','dlnR_fixedM')])
        replay = prod_pair(A,a,rho,rtol=1e-9,atol=1e-11,max_step=.05)
        grid.append(dict(alpha=a,A=A,rho_c=rho, M_std=std['M']/1.98847e33,
                         M_mod=mod['M']/1.98847e33,R_std_km=std['R']/1e5,
                         dlnR=float(independent[0]),dlnM=float(independent[1]),
                         dlnR_fixedM=float(independent[2]),fractional_fixedM=float(np.expm1(independent[2])),
                         saved_error=float(np.max(abs(independent-saved))),
                         saved_radius_error=float(abs(independent[0]-saved[0])),
                         saved_mass_error=float(abs(independent[1]-saved[1])),
                         saved_fixed_mass_error=float(abs(independent[2]-saved[2])),
                         replay_error=float(np.max(abs(replay-saved)))))
    checks['grid_has_108_rows'] = len(grid)==108
    checks['saved_grid_agrees_2e_6_absolute_log'] = max(r['saved_error'] for r in grid)<2e-6
    checks['production_replay_agrees_2e_8'] = max(r['replay_error'] for r in grid)<2e-8
    with (HERE/'grid_check.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(grid[0]));w.writeheader();w.writerows(grid)

    # Corners + fiducial + central boundary: tight integration, start radius,
    # production derivative cutoff, and explicit root at baseline target mass.
    cases=[(.0042,.02,1e5),(.0042,.02,3e6),(.0067,.5,1e5),(.0067,.5,3e6),
           (.0083,.15,1e6),(.0087,.15,1e6)]
    selected=[]
    for A,a,rho in cases:
        s,m=flux(0,a,rho),flux(A,a,rho)
        ref=shifts(s,m)
        tight=dict(rtol=2e-12,atol=2e-14,step=.003)
        conv=shifts(flux(0,a,rho,**tight),flux(A,a,rho,**tight))
        start=shifts(flux(0,a,rho,eps=1e-4),flux(A,a,rho,eps=1e-4))
        floor_results={}
        tight_floor_results={}
        original=production.THETA_FLOOR
        try:
            for floor in (1e-6,1e-8,1e-10):
                production.THETA_FLOOR=floor
                try:
                    floor_results[str(floor)]=prod_pair(A,a,rho).tolist()
                except RuntimeError as exc:
                    floor_results[str(floor)]={'failure':str(exc)}
                try:
                    tight_floor_results[str(floor)]=prod_pair(A,a,rho,rtol=2e-12,atol=2e-14,max_step=.003).tolist()
                except RuntimeError as exc:
                    tight_floor_results[str(floor)]={'failure':str(exc)}
        finally: production.THETA_FLOOR=original
        try:
            prod_tight=prod_pair(A,a,rho,rtol=2e-12,atol=2e-14,max_step=.003).tolist()
        except RuntimeError as exc:
            prod_tight={'failure':str(exc)}
        root=brentq(lambda lr:np.log(flux(A,a,np.exp(lr))['M']/s['M']),
                    np.log(rho)-.2,np.log(rho)+.2,xtol=1e-11)
        matched=flux(A,a,np.exp(root))
        # Exact standard scaling gives Rstd at Mmod, no linearization.
        rs_at_mmod=s['R']*(m['M']/s['M'])**(-1/3)
        identity=np.log(m['R']/rs_at_mmod)
        selected.append(dict(A=A,alpha=a,rho_c=rho,shifts=ref.tolist(),
            tighter_max_change=float(max(abs(conv-ref))),start_max_change=float(max(abs(start-ref))),
            production_floor_shifts=floor_results, production_tight=prod_tight,
            production_tight_floor_shifts=tight_floor_results,
            production_default_max_error=float(max(abs(np.array(floor_results['1e-08'])-ref))),
            exact_identity_error=float(abs(identity-ref[2])),
            target_mass_Msun=s['M']/1.98847e33, matched_rho=float(np.exp(root)),
            mass_root_relative_error=float(matched['M']/s['M']-1),
            fixed_baseline_mass_dlnR=float(np.log(matched['R']/s['R'])),
            mass_coordinate_change=float(np.log(matched['R']/s['R'])-ref[2])))
    checks['selected_convergence_2e_8'] = max(r['tighter_max_change'] for r in selected)<2e-8
    checks['selected_start_2e_8'] = max(r['start_max_change'] for r in selected)<2e-8
    checks['selected_default_production_2e_7'] = max(r['production_default_max_error'] for r in selected)<2e-7
    checks['fixed_mass_identity'] = max(r['exact_identity_error'] for r in selected)<1e-14
    checks['mass_root_accuracy'] = max(abs(r['mass_root_relative_error']) for r in selected)<1e-9
    files=[PROD/'model.py',PROD/'predict.py',PROD/'validate.py',PROD/'results_predictions.csv',Path(__file__)]
    result=dict(environment=dict(python=sys.version,numpy=np.__version__,scipy=scipy.__version__,platform=platform.platform()),
                inputs={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
                baseline=baseline,slopes=slopes,selected=selected,
                max_saved_grid_error=max(r['saved_error'] for r in grid),
                max_replay_error=max(r['replay_error'] for r in grid),
                grid_fractional_range=[min(r['fractional_fixedM'] for r in grid),max(r['fractional_fixedM'] for r in grid)],
                checks=checks,pass_all=all(checks.values()))
    (HERE/'results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ['max_saved_grid_error','max_replay_error','grid_fractional_range','checks','pass_all']},indent=2))
    if not result['pass_all']: raise SystemExit(1)

if __name__=='__main__':run()
