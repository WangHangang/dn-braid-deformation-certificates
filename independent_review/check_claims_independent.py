"""Independent randomized modular coefficient audit of supplied manuscript; no imports from source."""
import random,sys
P=1000000007
rng=random.Random(98765001)
def fmod(x):return x%P
def make(n,q):
 N=2*n; lam=pow(q,-1,P); d=(lam-q)%P
 rho={i:(n-i if i<=n else n+1-i) for i in range(1,N+1)}
 prime=lambda i:N+1-i
 B={(i,j):(q*(i==prime(j)) +d*(j>n)*(i==j) -d*(i>prime(j))*pow(q,-rho[i]-rho[j],P))%P for i in range(1,N+1) for j in range(1,N+1)}
 c={i:B[i,N] for i in range(1,N+1)}
 def Q(i,j):
  if i==j:return {(i,i):lam}
  if i==prime(j):return {(a,prime(a)):B[a][i] for a in []} # overridden below
  z={(j,i):1}
  if i>j:z[(i,j)]=d
  return z
 def qact(i,j):
  if i==prime(j):return {(a,prime(a)):B[(a,i)] for a in range(1,N+1) if B[(a,i)]}
  return Q(i,j)
 return N,lam,d,prime,B,c,qact

def trial(n,q):
 N,lam,d,prime,B,c,Q=make(n,q); I=range(2,N); V=list(range(1,N+1));
 # normalized H as if previous lemma: arbitrary zero weight block except first/N column restrictions;
 # s neg-neg=0, s_{n,N}=0, outer s_iN=0, H(ii)=0.
 s={};
 for i in V:
  for j in range(i+1,N+1):
   if j==prime(i):continue
   val=0 if (i>n and j>n) or j==N else rng.randrange(P)
   s[i,j]=val;s[j,i]=-val%P
 A={(i,j):rng.randrange(P) for i in V for j in V}
 for i in V:A[i,1]=0;A[i,N]=0
 gamma=rng.randrange(P);A[N,N]=gamma
 def H(i,j):
  if i==j:return {}
  if i==prime(j):return {(a,prime(a)):A[a,i] for a in V if A[a,i]}
  val=s[i,j]
  return {(j,i):val} if val else {}
 # apply right to left six monomials
 def apply(w,seq):
  vals={w:1}
  for op,ij in reversed(seq):
   nv={}
   for tup,x in vals.items():
    for (a,b),coef in (Q if op=='Q' else H)(tup[ij[0]],tup[ij[1]]).items():
     t=list(tup); t[ij[0]]=a; t[ij[1]]=b;t=tuple(t)
     nv[t]=(nv.get(t,0)+x*coef)%P
   vals={k:v for k,v in nv.items() if v}
  return vals
 p12=(0,1);p23=(1,2)
 seqs=[([('H',p12),('Q',p23),('Q',p12)],1),([('Q',p12),('H',p23),('Q',p12)],1),([('Q',p12),('Q',p23),('H',p12)],1),([('H',p23),('Q',p12),('Q',p23)],-1),([('Q',p23),('H',p12),('Q',p23)],-1),([('Q',p23),('Q',p12),('H',p23)],-1)]
 cache={}
 def L(inp,out):
  if inp not in cache:
   z={}
   for seq,sign in seqs:
    for t,v in apply(inp,seq).items():z[t]=(z.get(t,0)+sign*v)%P
   cache[inp]=z
  return cache[inp].get(out,0)
 errors=[];check=0
 def ck(name,inp,out,expected):
  nonlocal check
  check+=1
  a=L(inp,out);b=expected%P
  if a!=b:errors.append((name,inp,out,a,b))
 for i in I:
  ck('N11->1ii', (N,1,1),(1,i,prime(i)), -c[i]*(s[1,i]+s[1,prime(i)]))
  ck('N11->ii1', (N,1,1),(i,prime(i),1), d*sum(c[j]*A[i,j] for j in I)+lam*gamma*c[i])
 for k in I:
  kp=prime(k);b=s[1,k]
  ck('kN1->1kN',(k,N,1),(1,k,N),-c[kp]*A[1,k])
  ck('Nk1->Nk1',(N,k,1),(N,k,1),0) if False else None
  # output of rowN? manuscript unclear: Nk1 coefficient to (N,k,1)? stated for Nk1 terms 4-6 -ckp ck bk -ckp ANk -delta gamma; likely output k? search. 
  for i in I:
   if i in (k,kp):continue
   v=(s[k,i]+s[k,prime(i)])%P
   lhsA=(c[i]*(b+v)+d*c[kp]*(i<kp)*A[i,k])
   ck('A',(k,N,1),(i,prime(i),k),-lhsA)
   lhsC=(c[kp]*A[i,k]+d*c[i]*(i>kp)*v)
   ck('C',(k,N,1),(i,k,prime(i)),-lhsC)
  for r in V:
   term=(lam*c[k]*A[r,k]+c[kp]*B[(k,k)]*A[r,kp]+d*sum(c[j]*A[r,j] for j in I if j<k and j!=kp)+c[kp]*A[k,k]*B[(r,kp)]+b*c[k]*(r==kp))
   ck('E',(k,N,1),(k,r,prime(r)),-term)
 return check, errors
if __name__=='__main__':
 for n in (4,5):
  for q in (2,3):
   checks,errs=trial(n,q)
   print(f'n={n} q={q} coefficient checks={checks} mismatches={len(errs)}',flush=True)
   for row in errs[:12]:print(row,flush=True)
