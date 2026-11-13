"""Read-only audit of source-table arithmetic and cached score provenance.
Writes STEP09_INPUT_AUDIT.json only; does not import/run the scan or analyzer.
"""
from pathlib import Path
import ast, csv, hashlib, json, math
HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
MODEL=ROOT/'paper1_SolarCore/solarcore_v3/model.py'
ANALYZE=MODEL.with_name('analyze.py')
SCAN=MODEL.with_name('results_scan.csv')
def literal(path,name):
    for node in ast.parse(path.read_text()).body:
        if isinstance(node,ast.Assign) and any(isinstance(t,ast.Name) and t.id==name for t in node.targets):
            return ast.literal_eval(node.value)
    raise ValueError(name)
y=literal(MODEL,'EXP_R');sig=literal(MODEL,'EXP_SIGMA')
assert y==literal(ANALYZE,'EXP_R')
# Author manuscript Table 2, p.15: central, stat lower/upper, syst lower/upper,
# B16-GS98 Table 6 denominator and its fractional theory uncertainty.
inputs={'pp':(6.1,.5,.5,.5,.3,5.98,.006,10),
        'pep':(1.27,.19,.19,.12,.08,1.44,.01,8),
        'be7':(4.99,.11,.11,.08,.06,4.93,.06,9),
        'b8':(5.68,.41,.39,.03,.03,5.46,.12,6)}
rows=[]
for branch,(f,sl,su,el,eu,b,t,power) in inputs.items():
    lo=math.hypot(sl,el)/b;hi=math.hypot(su,eu)/b
    sym=math.sqrt((lo*lo+hi*hi)/2)
    rows.append(dict(branch=branch,flux_unit_power10_cm_minus2_s_minus1=power,
        measured_central=f,stat_lower=sl,stat_upper=su,syst_lower=el,syst_upper=eu,
        b16_gs98=b,theory_fraction=t,ratio=f/b,sigma_ratio_lower=lo,
        sigma_ratio_upper=hi,sigma_ratio_rms=sym,
        source_precision_diagonal_total=math.hypot(sym,t),legacy_ratio=y[branch],
        legacy_sigma_meas=sig[branch],legacy_sigma_total=math.hypot(sig[branch],t if branch in ('be7','b8') else 0)))
scan=list(csv.DictReader(SCAN.open()));residual=[]
for r in scan:
    recomputed=sum(((float(r['R_'+b])-y[b])/sig[b])**2 for b in y)
    residual.append(abs(recomputed-float(r['chi2'])))
null=[r for r in scan if float(r['A'])==0]
result={'sources':{'borexino_author_manuscript':'https://unige.iris.cineca.it/bitstream/11567/930098/5/11567_930098_merged_1534757340.pdf',
 'borexino_locator':'Table 2, manuscript page 15; visually checked. HZ pep and combined HE B8.',
 'b16':'https://arxiv.org/html/1611.09867v4','b16_locator':'Table 6; use GS98 model column, not Solar fit column.'},
 'source_rows':rows,'symmetrization':'sqrt((sigma_lower^2+sigma_upper^2)/2); chosen proxy approximation, not reconstructed likelihood',
 'legacy_null_measurement_score':sum(((1-y[b])/sig[b])**2 for b in y),
 'legacy_null_headline_score':sum(((1-r['legacy_ratio'])/r['legacy_sigma_total'])**2 for r in rows),
 'scan':{'rows':len(scan),'A_range':[min(float(r['A']) for r in scan),max(float(r['A']) for r in scan)],
 'alpha_range':[min(float(r['alpha']) for r in scan),max(float(r['alpha']) for r in scan)],
 'null_cached_score_values':sorted(set(float(r['chi2']) for r in null)),
 'cached_score_mismatch_count_at_1e_minus10':sum(v>1e-10 for v in residual),
 'max_cached_score_abs_difference':max(residual),
 'interpretation':'Cached chi2 is incompatible with current model.EXP_R/EXP_SIGMA; analyze.py recomputes scores from ratios.'},
 'evidence_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [MODEL,ANALYZE,SCAN,Path(__file__)]}}
(HERE/'STEP09_INPUT_AUDIT.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
