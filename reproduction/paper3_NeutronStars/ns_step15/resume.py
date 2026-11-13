"""Admit saved rows under unchanged scientific code; continue serially/throttled.
Original interrupted files remain immutable. A new manifest records adoption,
source equivalence, replay evidence, and the enforced resource settings.
"""
from pathlib import Path
import argparse,ast,json,platform,sys,traceback
import numpy as np
import scipy
from scan import ROOT,HERE,inputs,metrics,sha,SETTINGS,FRACTIONS,REFERENCES,EMBEDDINGS,summarize
from resource_control import require_controls,measurement

def admit(parent):
    require_controls();points,solar,ns=inputs();m=json.loads((parent/'manifest.json').read_text())
    if not (m['status']=='interrupted' or (m['status']=='failed' and m.get('error')=='KeyboardInterrupt()')):raise ValueError('Only explicitly interrupted source may be adopted')
    old=HERE/'diagnostics/scan_interrupted.py'
    for name,h in m['contract']['source_sha256'].items():
        p=ROOT/name
        if sha(p)!=h and Path(name).name in ['scan.py','resume.py','resource_control.py']:
            candidates=[HERE/'diagnostics/checkpoint25'/Path(name).name]
            if Path(name).name=='scan.py':candidates.append(old)
            p=next((q for q in candidates if q.exists() and sha(q)==h),p)
        if sha(p)!=h:raise ValueError('Interrupted source mismatch '+name)
    for p,h in m['contract']['solar_artifacts'].items():
        if sha(ROOT/p)!=h:raise ValueError('Solar source mismatch')
    for p,h in m['contract']['ns_audit_sources'].items():
        if sha(ROOT/p)!=h:raise ValueError('NS source mismatch')
    for k,v in [('python',platform.python_version()),('numpy',np.__version__),('scipy',scipy.__version__)]:
        if m['runtime'][k]!=v:raise ValueError('Runtime changed '+k)
    def science(path):
        tree=ast.parse(path.read_text());f=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='metrics')
        f.body=[n for n in f.body if not (isinstance(n,ast.ImportFrom) and n.module=='resource_control') and not (isinstance(n,ast.Expr) and isinstance(n.value,ast.Call) and isinstance(n.value.func,ast.Name) and n.value.func.id=='require_controls')]
        return ast.dump(f,include_attributes=False)
    if science(old)!=science(HERE/'scan.py'):raise ValueError('Scientific metrics implementation changed')
    if m['contract']['solver_settings']!=SETTINGS or m['contract']['points']!=points or m['contract']['fractions']!=FRACTIONS or m['contract']['reference_lengths_cm']!=REFERENCES or m['contract']['embeddings']!=EMBEDDINGS:raise ValueError('Scientific grid changed')
    jobs=[dict(p,A=p['boundary_A']*f,fraction=f,rstar_cm=ref,embedding=emb) for p in points for ref in REFERENCES for emb in EMBEDDINGS for f in FRACTIONS]
    ordered=[p for p in jobs if p['A']==0]+[p for p in jobs if p['A']>0]
    baseline=json.loads((parent/'baseline.json').read_text());gr=baseline['primary']
    rows=[]
    # Refuse even a malformed trailing line; never silently discard an ambiguous record.
    for line in (parent/'evaluations.jsonl').read_text().splitlines():rows.append(json.loads(line))
    if len(rows)!=m.get('saved_journal_lines',m.get('rows')):raise ValueError('Saved count mismatch')
    for row,point in zip(rows,ordered):
        if row['point']!=point:raise ValueError('Out-of-grid/duplicate/out-of-order saved row')
        if row['reused_GR']!=(point['A']==0):raise ValueError('Incorrect GR-reuse flag')
        if row['domain_peak'] not in row['evaluations'] or row['R14'] not in row['evaluations']:raise ValueError('Unbound output')
        if abs(row['R14']['M']-1.4)>5e-7 or len(row['pressure_grid'])!=13:raise ValueError('Invalid root/grid')
        for star in row['evaluations']:
            if not all(np.isfinite(v) and v>0 for v in star.values()):raise ValueError('Invalid stellar row')
            if star['P_surface']!=gr['R14']['P_surface']:raise ValueError('Surface mismatch')
        if any(v not in row['evaluations'] for v in row['pressure_grid']):raise ValueError('Unbound pressure grid')
        if row['domain_peak']['M']<max(v['M'] for v in row['pressure_grid'])-1e-9:raise ValueError('Missed saved maximum')
        if row['delta_mass']!=row['domain_peak']['M']/gr['domain_peak']['M']-1 or row['delta_R14']!=row['R14']['R']/gr['R14']['R']-1:raise ValueError('Shift algebra mismatch')
        if point['A']==0 and (row['domain_peak']!=gr['domain_peak'] or row['R14']!=gr['R14']):raise ValueError('Null row mismatch')
    if len(rows)>len(ordered):raise ValueError('Extra saved rows')
    return m,baseline,rows,ordered

def run(parent,out):
    require_controls();out.mkdir(parents=True,exist_ok=False)
    old,baseline,rows,jobs=admit(parent)
    files=[HERE/n for n in ['scan.py','resume.py','resource_control.py']]
    contract=dict(old['contract']);contract['source_sha256']=dict(old['contract']['source_sha256'])
    contract['source_sha256'].update({str(p.relative_to(ROOT)):sha(p) for p in files})
    m={'status':'running','contract':contract,'runtime':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,'workers':1},'adoption':{'parent_manifest_sha256':sha(parent/'manifest.json'),'parent_journal_sha256':sha(parent/'evaluations.jsonl'),'parent_baseline_sha256':sha(parent/'baseline.json'),'adopted_rows':len(rows),'scientific_metrics_AST_identical_excluding_resource_guard':True},'resource_controls':measurement()}
    def save():(out/'manifest.json').write_text(json.dumps(m,indent=2)+'\n')
    save();count=len(rows)
    try:
        # Replay common baseline and spread-out saved nonzero rows before accepting reuse.
        replay=[]
        for name,entry in baseline.items():
            fresh=metrics(entry['point'])
            if fresh!=entry:raise RuntimeError('Baseline replay differs '+name)
            replay.append({'baseline':name,'exact_match':True})
        ids=[i for i,r in enumerate(rows) if not r['reused_GR']]
        for i in [ids[0],ids[len(ids)//2],ids[-1]]:
            fresh=metrics(rows[i]['point']);oldrow=rows[i]
            for k in ['point','domain_peak','peak_status','R14','pressure_grid','evaluations']:
                if fresh[k]!=oldrow[k]:raise RuntimeError('Saved-row replay differs '+str(i)+' '+k)
            replay.append({'row_index':i,'point':fresh['point'],'exact_match':True})
        (out/'adoption_verification.json').write_text(json.dumps({'status':'pass','rows_structurally_verified':count,'replays':replay},indent=2)+'\n')
        (out/'baseline.json').write_text(json.dumps(baseline,indent=2)+'\n')
        gr=baseline['primary']
        with (out/'evaluations.jsonl').open('w') as f:
            for r in rows:f.write(json.dumps(r)+'\n')
            f.flush()
            print('Adopted',count,'verified rows;',len(jobs)-count,'remaining',flush=True)
            for p in jobs[count:]:
                value=metrics(p);value['delta_mass']=value['domain_peak']['M']/gr['domain_peak']['M']-1;value['delta_R14']=value['R14']['R']/gr['R14']['R']-1;value['reused_GR']=False
                f.write(json.dumps(value)+'\n');f.flush();rows.append(value)
                if len(rows)%10==0:print(len(rows),'of',len(jobs),'complete',flush=True)
        if len(rows)!=1980:raise RuntimeError('Incomplete grid')
        extrema={key:{q:rows.index(row) for q,row in val.items()} for key,val in summarize(rows).items()}
        (out/'summary.json').write_text(json.dumps({'rows':len(rows),'nonzero_rows':sum(not r['reused_GR'] for r in rows),'stellar_evaluations':sum(len(r['evaluations']) for r in rows if not r['reused_GR'])+sum(len(r['evaluations']) for r in baseline.values()),'extrema_indices':extrema,'global_envelope_certified':False},indent=2)+'\n')
        m.update(status='complete',rows=len(rows),resource_controls=measurement(),sha256={n:sha(out/n) for n in ['evaluations.jsonl','baseline.json','summary.json','adoption_verification.json']});save()
    except BaseException as exc:
        m.update(status='failed',error=repr(exc),traceback=traceback.format_exc(),rows=len(rows),resource_controls=measurement());save();raise
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();run(a.source,a.output)
