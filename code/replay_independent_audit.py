"""Replay supplied independent audit, and fill its omitted row-N diagnostic.

Public source files are never modified. Temporary copies are used for replay.
The public copies have had only an obsolete import path removed. SymPy is required for five symbolic
assertions; all modular checks themselves need only the standard library.
"""
from pathlib import Path
import contextlib,hashlib,io,json,re,runpy,subprocess,sys,tempfile


def main():
    root=Path(__file__).resolve().parents[1]
    supplied_dir=root/'independent_review'
    report={'supplied_script_directory':'anc/independent_review',
            'changes_to_received_originals':False,'source_scripts':{},'runs':[]}
    with tempfile.TemporaryDirectory(prefix='dn_review_replay_') as td:
        tmp=Path(td)
        for name in ['check_normal_form_symbolic.py','check_outer_independent.py','check_claims_independent.py']:
            raw=(supplied_dir/name).read_bytes();s=raw.decode('utf-8-sig')
            # Python already includes the executing script's own directory.
            portable='\n'.join(line for line in s.splitlines() if not line.startswith("sys.path.insert(0,'/mnt/data/"))+'\n'
            (tmp/name).write_text(portable,encoding='utf-8')
            report['source_scripts'][name]={'original_sha256':hashlib.sha256(raw).hexdigest(),
                                          'portability_change':'public scripts already omit obsolete sys.path insertion; calculations unchanged'}
        import sympy as S
        capture=io.StringIO()
        with contextlib.redirect_stdout(capture):ns=runpy.run_path(str(tmp/'check_normal_form_symbolic.py'))
        q=ns['q'];L=ns['L'];Q=ns['two_Q'];H=ns['two_H']
        p,u,v,d,xa,xb=[ns[k] for k in ['p','u','v','d','x_a','x_b']]
        checks=[(L((1,1,2),Q,H)[1,1,2],(q**-2-1)*p),
                (L((2,1,1),Q,H)[2,1,1],(1-q**-2)*((q*q+1)*xa-q*(u+v)-d)),
                (L((2,2,1),Q,H)[2,2,1],-(1-q**-2)*((q*q+1)*xb-q*(u+v)-d))]
        t=ns['tvars'];Q=ns['three_Q'];H=ns['three_H']
        checks.extend([(L((3,1,2),Q,H)[3,1,2],(q-q**-1)*(t[1,2]-t[1,3])),
                       (L((2,3,1),Q,H)[2,3,1],(q-q**-1)*(t[1,3]-t[2,3]))])
        for actual,expected in checks:assert S.cancel(actual-expected)==0
        report['symbolic_assertions']={'count':5,'sympy_version':S.__version__,'status':'PASS','stdout':capture.getvalue()}
        def run(name,kind):
            r=subprocess.run([sys.executable,str(tmp/name)],cwd=tmp,text=True,capture_output=True,check=True)
            print(r.stdout,end='',flush=True)
            if kind=='outer':
                rows=re.findall(r'^(\d+) (\d+) (\d+) errors \[\]$',r.stdout,re.M)
                assert len(rows)==6,r.stdout
                count=sum(int(row[2]) for row in rows)
            else:
                rows=re.findall(r'n=(\d+) q=(\d+) coefficient checks=(\d+) mismatches=(\d+)',r.stdout)
                assert len(rows)==4 and all(row[3]=='0' for row in rows),r.stdout
                count=sum(int(row[2]) for row in rows)
            report['runs'].append({'script':name,'scope':kind,'coefficient_checks':count,'status':'PASS','stdout':r.stdout})
            return count
        assert run('check_outer_independent.py','outer')==150
        assert run('check_claims_independent.py','received_internal')==628
        # The received internal script did not exercise Eq.(17), the row-N output.
        # Add exactly that coefficient; ignore other check names for this run.
        s=(tmp/'check_claims_independent.py').read_text(encoding='utf-8')
        needle="  ck('kN1->1kN',(k,N,1),(1,k,N),-c[kp]*A[1,k])"
        assert needle in s
        s=s.replace(needle,needle+"\n  ck('rowN',(k,N,1),(N,k,1),-(c[kp]*(A[N,k]+c[k]*b)+d*gamma))")
        s=s.replace(' def ck(name,inp,out,expected):\n  nonlocal check',
                    " def ck(name,inp,out,expected):\n  nonlocal check\n  if name!='rowN':return")
        (tmp/'check_rowN_added.py').write_text(s,encoding='utf-8')
        assert run('check_rowN_added.py','supplementary_row_N')==28
    report['received_coefficient_checks_replayed']=778
    report['additional_row_N_coefficient_checks']=28
    report['scope']='Independent finite-field/symbolic diagnostics; the uniform theorem is established by the written index proof, not these samples.'
    (root/'results/received_audit_replay.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print('RECEIVED AUDIT: 778 checks, 5 exact symbolic assertions, and 28 added row-N checks passed.',flush=True)

if __name__=='__main__':main()
