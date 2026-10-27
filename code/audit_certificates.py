"""Independently reconstruct supplied pivot minors. Standard library only.

The certificate merely chooses rows and columns: every Jacobian entry is rebuilt
from the matrix-unit Q formula and all six Frechet terms. No old verifier imports.
"""
from pathlib import Path
from collections import defaultdict
import argparse,hashlib,json,time
from check_uniform_identities import build_q


def evalmod(f,q,p,derivative=False):
    if not f:return 0
    return sum(c*(e[0] if derivative else 1)*pow(q,e[0],p) for e,c in f.terms)%p


def make_modq(n,q,p,derivative=False):
    return {src:{dst:z for dst,f in col.items() if (z:=evalmod(f,q,p,derivative))}
            for src,col in build_q(n).items()}


def modrows(src,Q,H,p):
    result={}
    words=[(1,[(Q,0),(Q,1),(H,0)]),(1,[(Q,0),(H,1),(Q,0)]),
           (1,[(H,0),(Q,1),(Q,0)]),(-1,[(Q,1),(Q,0),(H,1)]),
           (-1,[(Q,1),(H,0),(Q,1)]),(-1,[(H,1),(Q,0),(Q,1)])]
    for sign,word in words:
        state={(src,None):1}
        for op,pos in word:
            nxt={}
            for (t,label),c in state.items():
                for dst,v in op.get((t[pos],t[pos+1]),{}).items():
                    nt=t[:pos]+dst+t[pos+2:]
                    lab=v if op is H else label
                    z=c if op is H else c*v
                    key=nt,lab
                    nxt[key]=(nxt.get(key,0)+z)%p
            state={key:z for key,z in nxt.items() if z}
        for (t,label),c in state.items():
            row=result.setdefault(t,{})
            row[label]=(row.get(label,0)+sign*c)%p
    return result


def determinant(M,p):
    M=[row[:] for row in M];d=1;m=len(M)
    for j in range(m):
        k=next((i for i in range(j,m) if M[i][j]),None)
        if k is None:return 0
        if k!=j:M[k],M[j]=M[j],M[k];d=-d%p
        pv=M[j][j];d=d*pv%p;inv=pow(pv,-1,p)
        rowj=M[j]
        for i in range(j+1,m):
            if M[i][j]:
                row=M[i];factor=row[j]*inv%p
                for k in range(j+1,m):row[k]=(row[k]-factor*rowj[k])%p
                row[j]=0
    return d


def support(n):
    N=2*n;pairs=[(i,j) for i in range(1,N+1) for j in range(1,N+1)]
    def wt(i):return tuple(int(i==a)-int(N+1-i==a) for a in range(1,n+1))
    sums={src:tuple(x+y for x,y in zip(wt(src[0]),wt(src[1]))) for src in pairs}
    return sorted((src,dst) for src in pairs for dst in pairs if sums[src]==sums[dst])


def check_file(path):
    start=time.perf_counter();data=json.loads(path.read_text(encoding='utf-8'))
    n,p,q=data['n'],data['prime'],data['q_specialization'];Q=make_modq(n,q,p)
    S=support(n)
    assert S==[(tuple(a),tuple(b)) for a,b in data['support_indexing']]
    columns=data['independent_columns_indices']
    rows=[(tuple(a),tuple(b)) for a,b in data['pivot_rows_lexical']]
    assert len(rows)==len(columns)==data['linear_rank_mod_p']
    H={src:{} for src in Q}
    for j,k in enumerate(columns):
        src,dst=S[k];H[src][dst]=j
    inputs={src for src,dst in rows}
    rowvalues={src:modrows(src,Q,H,p) for src in inputs}
    matrix=[[rowvalues[src].get(dst,{}).get(j,0) for j in range(len(columns))] for src,dst in rows]
    h=hashlib.sha256()
    for row in matrix:
        for z in row:h.update(z.to_bytes(4,'little'))
    digest=h.hexdigest();det=determinant(matrix,p)
    assert det!=0
    if 'pivot_minor_determinant_mod_p' in data:
        assert det==data['pivot_minor_determinant_mod_p']
        assert digest==data['pivot_minor_rowmajor_sha256']
    result={'file':path.name,'input_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
            'n':n,'prime':p,'q':q,'ambient':len(S),'minor_size':len(rows),
            'independent_determinant':det,'rowmajor_sha256':digest,
            'matches_supplied_determinant':None if 'pivot_minor_determinant_mod_p' not in data else True,
            'elapsed_seconds':round(time.perf_counter()-start,3),'status':'PASS'}
    print(result,flush=True)
    return result


def tangent_check(n,q,p):
    """Verify correct logarithmic q derivative and expose old swapped entry."""
    Q=make_modq(n,q,p);D=make_modq(n,q,p,True)
    def derivative_residual(D):
        # Treat its coefficients as labels, then sum; modrows expects a label,
        # so each transition is assigned an integer index with coefficient list.
        H={src:{} for src in Q};coeff=[]
        for src,col in D.items():
            for dst,z in col.items():H[src][dst]=len(coeff);coeff.append(z)
        N=2*n
        for a in range(1,N+1):
            for b in range(1,N+1):
                for c in range(1,N+1):
                    for dst,row in modrows((a,b,c),Q,H,p).items():
                        val=sum(z*coeff[j] for j,z in row.items())%p
                        if val:return {'input':[a,b,c],'output':list(dst),'value':val}
        return None
    good=derivative_residual(D);assert good is None
    # Old script used swap target (i,j), rather than unchanged (j,i).
    bad={src:col.copy() for src,col in D.items()};N=2*n
    v=(-pow(q,-1,p)-q)%p
    for i in range(1,N+1):
        for j in range(i+1,N+1):
            for dst,z in [((j,i),-v),((i,j),v)]:
                bad[j,i][dst]=(bad[j,i].get(dst,0)+z)%p
    fail=derivative_residual(bad);assert fail is not None
    return {'n':n,'q':q,'prime':p,'correct_q_derivative_residual':None,
            'old_q_derivative_first_nonzero':fail,'status':'PASS'}


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('certificate_directory',type=Path)
    ap.add_argument('--output',type=Path,default=Path(__file__).resolve().parents[1]/'results'/'independent_certificate_audit.json')
    args=ap.parse_args()
    results=[check_file(path) for path in sorted(args.certificate_directory.glob('*.json'))]
    tangents=[tangent_check(4,2,1000000007)]
    report={'method':'independent matrix-unit construction, six-term Jacobian, exact modular determinant',
            'certificates':results,'tangent_diagnostics':tangents}
    args.output.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print('all certificate and derivative checks passed',flush=True)

if __name__=='__main__':main()
