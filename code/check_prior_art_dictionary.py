"""Exact coefficient dictionary with Isaev--Ogievetsky--Pyatov (2004), Eq.43.

All q, alternating-twist, basis and scalar parameters are independent Laurent
variables. Finite-rank diagnostics accompany the uniform reindexing in the note.
Python standard library only; no floating point or numerical sampling.
"""
from itertools import combinations
from pathlib import Path
import json
from laurent_core import variable,accum


def check(n):
    N=2*n;ip=lambda i:N+1-i;rho=lambda i:n-i if i<=n else n+1-i
    pairs=list(combinations(range(n+1),2));nv=(n+1)*(n+2)//2
    q=variable(nv,0);one=q**0;delta=q**-1-q
    twist={uv:variable(nv,k+1) for k,uv in enumerate(pairs)}
    diagonals={i:variable(nv,1+len(pairs)+i-1) for i in range(1,n)}
    scalar=variable(nv,nv-1)
    def mu(i):return [int(i==a)-int(ip(i)==a) for a in range(1,n+1)]+[1]
    def chi(i,j):
        x,y=mu(i),mu(j);z=one
        for (u,v),t in twist.items():z=z*t**(x[u]*y[v]-x[v]*y[u])
        return z
    def D(i):return diagonals.get(i,one)
    def eta(i):
        x=mu(i);z=one
        for a in range(n):z=z*twist[a,n]**x[a]
        return z
    # Uniform dictionary: s=q^-1, reflected basis J, f_ij=chi(mu_i',mu_j').
    def f(i,j):return chi(ip(i),ip(j))
    condition_checks=0
    for i in range(1,N+1):
        assert not f(ip(i),i)*f(i,ip(i))**-1-eta(i)**4
        condition_checks+=1
        for j in range(1,N+1):
            assert not f(i,j)*f(ip(i),j)-eta(j)**2
            assert not f(i,j)*f(i,ip(j))-eta(i)**-2
            condition_checks+=2
    # Our operator is reconstructed from its weight-block actions.
    ours={}
    def add_ours(src,dst,z):
        ratio=chi(*dst)*chi(*src)**-1
        ratio=ratio*D(dst[0])*D(dst[1])*(D(src[0])*D(src[1]))**-1
        accum(ours,(src,dst),scalar*z*ratio)
    for a in range(1,N+1):
        for b in range(1,N+1):
            if a==b:add_ours((a,b),(a,b),q**-1)
            elif b!=ip(a):
                add_ours((a,b),(b,a),one)
                if a>b:add_ours((a,b),(a,b),delta)
            else:
                for i in range(1,N+1):
                    z=q*int(i==ip(a))+delta*int(a>n and i==a)
                    if i>ip(a):z=z-delta*q**(-rho(i)-rho(a))
                    add_ours((a,b),(i,ip(i)),z)
    # Literal matrix-unit formula Eq.(43), not the code constructing ours.
    lit={};s=q**-1;ell=s-s**-1
    def E(a,b,c,d,z):accum(lit,((b,d),(a,c)),z)
    for i in range(1,N+1):
        for j in range(1,N+1):
            E(i,j,j,i,s**(int(i==j)-int(i==ip(j)))*f(i,j)*f(j,i)**-1)
            if i<j:E(i,i,j,j,ell)
            if i>j:E(ip(i),j,i,ip(j),-ell*s**(rho(i)-rho(j))*f(ip(i),i)*f(j,ip(j))**-1)
    transformed={}
    for (src,dst),z in lit.items():
        a,b=tuple(map(ip,src)),tuple(map(ip,dst))
        ratio=D(b[0])*D(b[1])*(D(a[0])*D(a[1]))**-1
        accum(transformed,(a,b),scalar*ratio*z)
    for key in set(ours)|set(transformed):
        assert not ours.get(key,0)-transformed.get(key,0),(n,key)
    return {'n':n,'independent_Laurent_parameters':nv,
            'coefficient_positions_compared':len(set(ours)|set(transformed)),
            'twisting_constraint_identities':condition_checks,
            'result':'PASS'}


if __name__=='__main__':
    results=[check(n) for n in (4,5)]
    report={'source':'Isaev--Ogievetsky--Pyatov, Proc. Steklov Inst. Math.246 (2004),134-141; arXiv:math/0509251, Eq.(43)',
            'dictionary':'s=q^-1; J(v_i)=v_i_prime; f_ij=chi(mu_i_prime,mu_j_prime); external reflected diagonal basis; same scalar',
            'scope':'Exact finite-rank diagnostics of the written uniform change of indices, not a proof of novelty or all-rank exhaustion',
            'results':results}
    (Path(__file__).resolve().parents[1]/'results/prior_art_dictionary_checks.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report,indent=2))
