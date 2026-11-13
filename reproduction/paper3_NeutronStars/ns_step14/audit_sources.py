"""Verify pinned local EOS data against supplied official archive and source PDF.
Download URLs are recorded in output; no remote file is silently substituted.
Requires pdftotext. Inputs remain separate from the scientific reproduction.
"""
from pathlib import Path
import argparse,hashlib,json,re,subprocess,zipfile
import numpy as np
from model import ROOT,EOS,C

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def run(archive,pdf,output):
    z=zipfile.ZipFile(archive);matches={}
    for name in ['eos.nb','eos.thermo','eos.mr']:
        local=ROOT/'paper3_NeutronStars/ns_v7'/name
        matches[name]={'sha256':sha(local),'identical_to_official_archive':z.read(name)==local.read_bytes()}
        if not matches[name]['identical_to_official_archive']:raise ValueError('source mismatch '+name)
    text=subprocess.check_output(['pdftotext','-layout',str(pdf),'-'],text=True)
    segments=[text.split('Table 3. Equation of state of the inner crust.')[1].split('is 1.35')[0],text.split('Table 5. Equation of state of the liquid neutron-star core')[1].split('However,')[0]]
    rows=[]
    for seg in segments:
        for line in seg.splitlines():
            clean=re.sub(r'\s+E\s*([+-]?\d+)',r'E\1',line)
            try:values=list(map(float,clean.split()))
            except ValueError:continue
            if len(values)==8:
                rows.extend([values[:3],values[4:7]])
    if len(rows)!=78:raise ValueError('unexpected source table parse')
    a=np.loadtxt(ROOT/'paper3_NeutronStars/037/eos_sly4_dh2001_v2.dat')
    matched=[]
    for rho,p,eps,nb in a:
        hits=[v for v in rows if np.allclose([nb/1e39,rho,p],v,rtol=1e-10,atol=0)]
        if len(hits)!=1:raise ValueError('unmatched or ambiguous DH row '+str(rho))
        matched.append(hits[0])
    unused=[v for v in rows if v not in matched]
    # Compare independent loader with historical loader at all knots, including units.
    import sys
    sys.path.insert(0,str(ROOT/'paper3_NeutronStars/ns_v7'))
    from eos import CompOSE_SLy4
    fresh=EOS();old=CompOSE_SLy4()
    diffs={k:float(np.max(np.abs(getattr(fresh,k)/getattr(old,k)-1))) for k in ['P','eps','rho']}
    if max(diffs.values())>1e-14:raise ValueError('loader mismatch')
    result={'status':'pass','official_archive_url':'https://compose.obspm.fr/download/1D/NS/Skyrme/SLY4/eos.zip','official_archive_sha256':sha(archive),'official_members':matches,'DH_PDF_url':'https://arxiv.org/pdf/astro-ph/0111092','DH_PDF_sha256':sha(pdf),'DH_source_rows':len(rows),'DH_local_rows_matched':len(matched),'DH_omitted_source_rows_nb_rho_P':unused,'DH_wrong_epsilon_over_rho_c2_range':[float(np.min(a[:,2]/(a[:,0]*C*C))),float(np.max(a[:,2]/(a[:,0]*C*C)))],'independent_vs_legacy_compose_loader_max_relative_differences':diffs,'compose_manual':'https://compose.obspm.fr/download/pdf/manual_v2.00.pdf','compose_formulas':{'pressure':'Q1 * nb * MeV_fm3_to_cgs','energy_density':'mn * (1+Q7) * nb * MeV_fm3_to_cgs','law_density':'mu * nb * 1e39; mu=1.66053906660e-24 g'},'source_sha256':{str(Path(__file__).resolve().relative_to(ROOT)):sha(Path(__file__))}}
    output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--compose-zip',type=Path,required=True);p.add_argument('--dh-pdf',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();run(a.compose_zip,a.dh_pdf,a.output)
