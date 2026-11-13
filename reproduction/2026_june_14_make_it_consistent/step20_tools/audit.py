"""Small independent algebra/arithmetic controls; no stellar integrations."""
from pathlib import Path
import ast,hashlib,json,math,os,subprocess
from resource_control import require_controls,measurement
require_controls()
ROOT=Path(__file__).resolve().parents[2];E=ROOT/'2026_june_14_make_it_consistent/step20_evidence'
# Live inventory includes this audit and the one short-lived ps inspector, no pool.
ps=subprocess.run(['/bin/ps','-axo','pid,ppid,ni,%cpu,command'],capture_output=True,text=True,check=True).stdout
pid=os.getpid();lines=ps.splitlines();own=[]
for line in lines[1:]:
 fields=line.split(None,4)
 if len(fields)==5 and (int(fields[0])==pid or int(fields[1])==pid):own.append(line)
assert any(int(l.split()[0])==pid and int(l.split()[2])==19 for l in own)
assert all(int(l.split()[0])==pid or '/bin/ps' in l for l in own)
(E/'live_process.txt').write_text(lines[0]+'\n'+'\n'.join(own)+'\n')
A=.0087
checks=[]
for alpha in (.02,.12,.334,.5,1):
 for density in (.01,1,100,1e8):
  x=density**alpha;s=A*x/(1+x);up=(s/(A-s))**(1/alpha)
  assert math.isclose(up,density,rel_tol=16*math.ulp(1.)/(alpha*(1-s/A)))
  # Algebraic inversion residual vs actual weak-curvature scalar residual.
  checks.append({'alpha':alpha,'rho_over_rho0':density,'sigma':s,'assigned_equilibrium_residual':up-density,'action_source_residual_over_rho0':up+density/2})
  assert checks[-1]['action_source_residual_over_rho0']>0
# Exact rational trace elimination, arbitrary positive constants and potential.
from fractions import Fraction as Q
mp2=Q(7);kin=Q(3);s=Q(1,100);F=mp2*(1-s);T=-Q(5);U=Q(2);up=Q(4);grad2=Q(1,9)
box=(up+mp2*(-T+kin*grad2+4*U)/(2*F))/(kin+3*mp2**2/(2*F))
R=(-T+kin*grad2+4*U-3*mp2*box)/F
assert kin*box==up+mp2*R/2
# Regular nonuniform density, constant deficit: analytic weighted integral.
d=-.02;q=.7;t=1
op=(1/(3+d)-q*t*t/(5+d))/(1/3-q*t*t/5);local=3/(3+d)
assert abs(op-local)>1e-4
# Geometry: Schwarzschild shell slope is nonzero despite zero density.
geom=3-1/(10-2);assert geom==2.875
# NIST covariance propagation for descriptive fixed-reference contrasts.
data=json.loads((E/'nist_tables.json').read_text());g=data['G_in_1e_minus_11_SI'];u=[g[i]*data['relative_standard_uncertainty_ppm_table18'][i]*1e-6 for i in range(4)]
cov=[[u[i]*u[j]*data['correlation'][i][j] for j in range(4)] for i in range(4)]
# Cholesky positivity without a numerical-library dependency.
L=[[0.]*4 for _ in range(4)]
for i in range(4):
 for j in range(i+1):
  val=cov[i][j]-sum(L[i][k]*L[j][k] for k in range(j))
  L[i][j]=math.sqrt(val) if i==j else val/L[j][j]
contrasts=[]
for name,c in [('copper_minus_sapphire_servo',[1,0,-1,0]),('copper_minus_sapphire_free',[0,1,0,-1]),('free_minus_servo_copper',[-1,1,0,0]),('free_minus_servo_sapphire',[0,0,-1,1])]:
 val=sum(c[i]*g[i] for i in range(4));var=sum(c[i]*cov[i][j]*c[j] for i in range(4) for j in range(4))
 contrasts.append({'name':name,'ppm_at_fixed_Gref_6_6737':val/6.6737*1e6,'standard_uncertainty_ppm_excludes_dark_bias':math.sqrt(var)/6.6737*1e6})
# Exact vs linear material ansatz. These are chosen inputs, not measured densities.
pairs=[]
for alpha in (.02,.12,.334,.5):
 def sig(r):
  x=(r/150)**alpha;return A*x/(1+x)
 s1,s2=sig(19.3),sig(2.7)
 pairs.append({'alpha':alpha,'leading_ppm':(s1-s2)*1e6,'exact_ppm':((1-s2)/(1-s1)-1)*1e6})
# Confirm old script's actual parameter differs from its manuscript's .12.
tree=ast.parse((ROOT/'paper6/geff_check.py').read_text());found=False
for node in tree.body:
 if isinstance(node,ast.Assign) and isinstance(node.targets[0],ast.Tuple):
  if [n.id for n in node.targets[0].elts]==['A','ALPHA']:
   assert ast.literal_eval(node.value)==(.0087,.122);found=True
assert found
out={'status':'pass','new_stellar_runs':0,'source_sign_checks':checks,'trace_elimination_exact':True,'nonlocal_counterexample':{'operator_ratio':op,'local_central_ratio':local},'schwarzschild_dimension_at_r_over_GM_10':geom,'nist_contrasts':contrasts,'nist_covariance_positive_definite':True,'material_ansatz':pairs,'historical_local_alpha':.122,'resource':measurement()}
(E/'audit.json').write_text(json.dumps(out,indent=2)+'\n')
# Generated tables are the only new numerical manuscript assets.
def table(headers,rows):
 return '\\begin{tabular}{'+'l'*len(headers)+'}\n\\toprule\n'+' & '.join(headers)+' \\\\\n\\midrule\n'+''.join(' & '.join(row)+' \\\\\n' for row in rows)+'\\bottomrule\n\\end{tabular}\n'
(E/'material_table.tex').write_text(table([r'$\alpha$',r'Leading (ppm)',r'Exact (ppm)'],[[str(p['alpha']),f"{p['leading_ppm']:.2f}",f"{p['exact_ppm']:.2f}"] for p in pairs]))
labels=['Cu--Sa, servo','Cu--Sa, free','Free--servo, Cu','Free--servo, Sa']
(E/'contrast_table.tex').write_text(table(['Contrast','Difference (ppm)',r'$u$ (ppm)'],[[labels[i],f"{p['ppm_at_fixed_Gref_6_6737']:.2f}",f"{p['standard_uncertainty_ppm_excludes_dark_bias']:.2f}"] for i,p in enumerate(contrasts)]))
for f in ('material_table.tex','contrast_table.tex'):(ROOT/'paper6/terrestrial_G_note'/f).write_bytes((E/f).read_bytes())
print(json.dumps({'controls':'pass','contrasts':contrasts,'operator_counterexample':out['nonlocal_counterexample']},indent=2))
