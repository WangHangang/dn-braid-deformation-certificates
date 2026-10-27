import sys,random
from check_claims_independent import make,P
R=random.Random(12081211)
def test(n,q):
 N,lam,d,ip,B,c,Q=make(n,q);I=range(2,N);V=range(1,N+1)
 s={}
 for i in V:
  for j in range(i+1,N+1):
   if j==ip(i):continue
   val=0 if (i>n and j>n) or (j==N and i==n) else R.randrange(P)
   s[i,j]=val;s[j,i]=-val%P
 A={(i,j):R.randrange(P) for i in V for j in V}
 for i in I:A[i,N]=0
 def H(a,b):
  if a==b:return {}
  if a==ip(b):return {(i,ip(i)):A[i,a] for i in V if A[i,a]}
  return {(b,a):s[a,b]} if s[a,b] else {}
 p12=(0,1);p23=(1,2)
 seqs=[([('H',p12),('Q',p23),('Q',p12)],1),([('Q',p12),('H',p23),('Q',p12)],1),([('Q',p12),('Q',p23),('H',p12)],1),([('H',p23),('Q',p12),('Q',p23)],-1),([('Q',p23),('H',p12),('Q',p23)],-1),([('Q',p23),('Q',p12),('H',p23)],-1)]
 def L(inp,out):
  answer=0
  for seq,sgn in seqs:
   states={inp:1}
   for op,ij in reversed(seq):
    nxt={}
    for w,z in states.items():
     for pair,fac in (Q if op=='Q' else H)(w[ij[0]],w[ij[1]]).items():
      t=list(w);t[ij[0]],t[ij[1]]=pair;t=tuple(t)
      nxt[t]=(nxt.get(t,0)+z*fac)%P
    states=nxt
   answer=(answer+sgn*states.get(out,0))%P
  return answer
 err=[]; count=0
 def ck(stage,src,dst,exp):
  nonlocal count
  count+=1
  val=L(src,dst)
  if val!=exp%P:err.append((stage,src,dst,val,exp%P))
 ck('firstdiag',(1,1,N),(1,1,N),(lam**2-q*q)*A[1,1])
 for i in I:ck('firstcol',(1,1,N),(i,1,ip(i)),d*A[i,1])
 # post forced col1: 
 A[1,1]=0
 for i in I:A[i,1]=0
 for i in I:ck('second',(N,N,1),(i,ip(i),N),c[i]*(lam*A[1,N]+s[i,N]+s[ip(i),N]))
 # post forced zeta & s_iN
 A[1,N]=0
 for i in I:s[i,N]=0;s[N,i]=0
 for i in I:ck('third',(N,1,N),(N,i,ip(i)),-lam*c[i]*A[N,1])
 return count,err
for n in (4,5,6):
 for q in (2,3):
  count,errors=test(n,q)
  print(n,q,count,'errors',errors[:3])
