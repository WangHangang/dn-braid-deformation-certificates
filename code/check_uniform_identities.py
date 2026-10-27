"""Exact Laurent checks of the coefficient identities used in the uniform proof.

These finite-rank checks are diagnostics, not a replacement for the index-uniform
coefficient expansion in the proof. Python standard library only.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
from time import perf_counter
from laurent_core import Laurent, variable, accum

q = variable(1, 0)
one = q ** 0
lam = q ** -1
delta = lam - q


def rho(n, i):
    return n - i if i <= n else n + 1 - i


def build_q(n):
    """Independent transcription of matrix-unit formula for Q=PR (1-based)."""
    N = 2*n
    out = {(i,j): {} for i in range(1,N+1) for j in range(1,N+1)}
    def add(src,dst,z):
        accum(out[src],dst,z)
    for i in range(1,N+1):
        add((i,i),(i,i),lam)
        add((i,N+1-i),(N+1-i,i),q)
        for j in range(1,N+1):
            if i != j and i != N+1-j:
                add((i,j),(j,i),one)
            if i < j:
                add((j,i),(j,i),delta)
            if i > j:
                add((N+1-j,j),(i,N+1-i),-delta*q**(rho(n,j)-rho(n,i)))
    return out


def lfadd(*forms):
    out = {}
    for f in forms:
        for k,v in f.items():
            accum(out,k,v)
    return out


def scale(c,f):
    return {k:c*v for k,v in f.items() if c*v}


def h_builder(n,stage):
    N=2*n
    I=range(2,N)
    def S(a,b):
        if a==b or a+b==N+1 or (a>n and b>n): return {}
        if {a,b}=={n,N}: return {}
        if stage >= 1 and N in (a,b): return {}
        return {('s',min(a,b),max(a,b)):one if a<b else -one}
    def A(i,j):
        if stage==0:
            if j==N and i in I:return {}
        elif stage==1:
            if j==1 and i!=N:return {}
            if j==N and i!=N:return {}
        else:
            if j==1:return {}
            if j==N and i!=N:return {}
        return {('A',i,j):one}
    out={(a,b):{} for a in range(1,N+1) for b in range(1,N+1)}
    for a,b in out:
        if b==N+1-a:
            for i in range(1,N+1):
                f=A(i,a)
                if f:out[a,b][i,N+1-i]=f
        else:
            f=S(a,b)
            if f:out[a,b][b,a]=f
    return out,S,A


def linearized_rows(src,Q,H):
    # State coefficients carry either no H-label, or exactly one label.
    result={}
    words=[(1,[(Q,0),(Q,1),(H,0)]),
           (1,[(Q,0),(H,1),(Q,0)]),
           (1,[(H,0),(Q,1),(Q,0)]),
           (-1,[(Q,1),(Q,0),(H,1)]),
           (-1,[(Q,1),(H,0),(Q,1)]),
           (-1,[(H,1),(Q,0),(Q,1)])]
    for sign,word in words:
        state={(src,None):one}
        for op,pos in word:
            nxt={}
            for (t,label),c in state.items():
                for dst,v in op[t[pos],t[pos+1]].items():
                    nt=t[:pos]+dst+t[pos+2:]
                    if op is H:
                        assert label is None
                        for lab,co in v.items():accum(nxt,(nt,lab),c*co)
                    else:accum(nxt,(nt,label),c*v)
            state=nxt
        for (t,label),c in state.items():
            if t not in result:result[t]={}
            accum(result[t],label,sign*c)
    return result


def check_rank(n):
    N=2*n; I=range(2,N); ip=lambda k:N+1-k
    Q=build_q(n)
    B=lambda i,j:Q[j,ip(j)].get((i,ip(i)),0)
    c=lambda i:B(i,N)
    count=0
    def check(name,got,want):
        nonlocal count
        dif=lfadd(got,scale(-one,want))
        if dif:
            serial={str(k):v.serial() for k,v in dif.items()}
            raise AssertionError((n,name,serial))
        count+=1
    # Direct zero-block formula, including all exceptional corners.
    for i in range(1,N+1):
        for j in range(1,N+1):
            v=q*int(i==ip(j))+delta*int(j>n and i==j)
            if i>ip(j):v=v-delta*q**(-rho(n,i)-rho(n,j))
            assert not B(i,j)-v, (n,i,j)
            count+=1
    H,S,A=h_builder(n,0)
    rows=linearized_rows((1,1,N),Q,H)
    check('outer_col1_diag',rows.get((1,1,N),{}),scale(lam**2-q**2,A(1,1)))
    for i in I:
        check(('outer_col1',i),rows.get((i,1,ip(i)),{}),scale(delta,A(i,1)))
    # First use the forced zeros in column 1 before examining (N,N,1).
    for dst in H[1,N].copy():
        if dst!=(N,1): del H[1,N][dst]
    rows=linearized_rows((N,N,1),Q,H)
    for i in I:
        want=scale(c(i),lfadd(scale(lam,A(1,N)),S(i,N),S(ip(i),N)))
        check(('outer_swap',i),rows.get((i,ip(i),N),{}),want)
    H,S,A=h_builder(n,1)
    rows=linearized_rows((N,1,N),Q,H)
    for i in I:
        check(('outer_eta',i),rows.get((N,i,ip(i)),{}),scale(-lam*c(i),A(N,1)))
    H,S,A=h_builder(n,2)
    rows=linearized_rows((N,1,1),Q,H)
    for i in I:
        check(('b_pair',i),rows.get((1,i,ip(i)),{}),scale(-c(i),lfadd(S(1,i),S(1,ip(i)))))
        want=lfadd(*(scale(delta*c(j),A(i,j)) for j in I),scale(lam*c(i),A(N,N)))
        check(('D',i),rows.get((i,ip(i),1),{}),want)
    for k in I:
        rows=linearized_rows((k,N,1),Q,H)
        for i in I:
            if i==k or i==ip(k):continue
            want=lfadd(scale(c(i),lfadd(S(1,k),S(k,i),S(k,ip(i)))),
                       scale(delta*c(ip(k))*int(i<ip(k)),A(i,k)))
            check(('A',k,i),rows.get((i,ip(i),k),{}),scale(-one,want))
            want=lfadd(scale(c(ip(k)),A(i,k)),
                       scale(delta*c(i)*int(i>ip(k)),lfadd(S(k,i),S(k,ip(i)))))
            check(('B',k,i),rows.get((i,k,ip(i)),{}),scale(-one,want))
        check(('row1',k),rows.get((1,k,N),{}),scale(-c(ip(k)),A(1,k)))
        want=lfadd(scale(c(ip(k)),lfadd(A(N,k),scale(c(k),S(1,k)))),scale(delta,A(N,N)))
        check(('rowN',k),rows.get((N,k,1),{}),scale(-one,want))
        for r in range(1,N+1):
            want=lfadd(scale(lam*c(k),A(r,k)),
                       scale(c(ip(k))*B(k,k),A(r,ip(k))),
                       *(scale(delta*c(j),A(r,j)) for j in I if j<k and j!=ip(k)),
                       scale(c(ip(k))*B(r,ip(k)),A(k,k)),
                       scale(c(k)*int(r==ip(k)),S(1,k)))
            check(('E',k,r),rows.get((k,r,ip(r)),{}),scale(-one,want))
    return {'n':n,'coefficient_identities_checked':count,'result':'PASS'}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--ranks',type=int,nargs='+',default=[4,5,6])
    parser.add_argument('--output',type=Path,default=Path(__file__).resolve().parents[1]/'results'/'uniform_identity_checks.json')
    args=parser.parse_args();start=perf_counter()
    results=[]
    for n in args.ranks:
        result=check_rank(n);results.append(result);print(result,flush=True)
    report={'method':'exact Laurent expansion in Z[q,q^-1]; no sampling',
            'scope':'diagnostic checks of the written uniform coefficient proof; finite checks alone do not prove the theorem',
            'results':results,'elapsed_seconds':perf_counter()-start}
    args.output.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')

if __name__=='__main__':main()
