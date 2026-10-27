"""Recheck universal two/three-letter normal forms, without SymPy."""
from pathlib import Path
import json
from itertools import combinations,permutations,product
from check_uniform_identities import build_q,linearized_rows,lfadd,scale,q,one,lam


def unit(label):return {label:one}


def main():
    Q=build_q(4)
    H={k:{} for k in Q}
    a,b,c=2,3,4
    H[a,a][a,a]=unit('xa');H[b,b][b,b]=unit('xb')
    H[a,b][a,b]=unit('p');H[a,b][b,a]=unit('u')
    H[b,a][a,b]=unit('v');H[b,a][b,a]=unit('d')
    report={}
    def ser(f):return {str(k):v.serial() for k,v in f.items()}
    for src in [(a,a,b),(b,a,a),(b,b,a)]:
        f=linearized_rows(src,Q,H).get(src,{})
        report[str(src)]=ser(f)
    # The first equation kills p; then the next two give the diagonal formulas.
    assert linearized_rows((a,a,b),Q,H)[a,a,b]=={'p':lam**2-one}
    H[a,b].pop((a,b))
    for src,x,sgn in [((b,a,a),'xa',1),((b,b,a),'xb',-1)]:
        expected=scale(sgn*(one-lam**2),lfadd(scale(q*q+one,unit(x)),scale(-q,lfadd(unit('u'),unit('v'))),scale(-one,unit('d'))))
        assert not lfadd(linearized_rows(src,Q,H)[src],scale(-one,expected))
    H={k:{} for k in Q}
    for i in (a,b,c):H[i,i][i,i]=unit('x')
    for i,j in combinations((a,b,c),2):
        H[i,j][j,i]=unit(('u',i,j))
        H[j,i][i,j]=unit(('v',i,j))
        H[j,i][j,i]=lfadd(scale(q*q+one,unit('x')),scale(-q,lfadd(unit(('u',i,j)),unit(('v',i,j)))))
    for src in [(c,a,b),(b,c,a)]:
        f=linearized_rows(src,Q,H).get(src,{})
        report[str(src)]=ser(f)
        print(src,ser(f))
        left,right=((a,b),(a,c)) if src==(c,a,b) else ((a,c),(b,c))
        difference=lfadd(unit(('u',*left)),unit(('v',*left)),
                         scale(-one,lfadd(unit(('u',*right)),unit(('v',*right)))))
        assert not lfadd(f,scale(-(q-lam),difference))
    # Check the entire three-letter system is zero after imposing one common t.
    for i,j in combinations((a,b,c),2):
        H[j,i][i,j]=lfadd(unit('t'),scale(-one,unit(('u',i,j))))
        H[j,i][j,i]=lfadd(scale(q*q+one,unit('x')),scale(-q,unit('t')))
    for src in permutations((a,b,c)):
        assert all(not f for f in linearized_rows(src,Q,H).values())
    report['status']='PASS: correct two-letter formulas and common exchange sum suffice on all three-letter permutations'
    path=Path(__file__).resolve().parents[1]/'results'/'normal_form_checks.json'
    path.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(report['status'])

if __name__=='__main__':main()
