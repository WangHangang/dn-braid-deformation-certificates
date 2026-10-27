import sympy as S
q=S.symbols('q',nonzero=True);l=1/q;delta=l-q
p,u,v,d,x_a,x_b=S.symbols('p u v d x_a x_b')
def two_Q(a,b):
 if a==b:return {(a,b):l}
 return {(b,a):S.Integer(1),**({(a,b):delta} if a>b else {})}
def two_H(a,b):
 if a==b:return {(a,b):(x_a if a==1 else x_b)}
 return {(1,2):p,(2,1):u} if a<b else {(1,2):v,(2,1):d}
def L(inp,Q,H):
 opseqs=[('HQQ',1,'121'),('QHQ',1,'121'),('QQH',1,'121'),('HQQ',-1,'232'),('QHQ',-1,'232'),('QQH',-1,'232')]
 result={}
 for word,sgn,loc in opseqs:
  states={inp:1}
  for op,place in reversed(list(zip(word,loc))):
   pos=0 if place=='1' or place=='2' and loc=='121' else 1
   # intentional explicit define using paired location sequences instead below
  ops=([('H',0),('Q',1),('Q',0)] if word=='HQQ' and loc=='121' else None)
  if loc=='121':slots=[0,1,0]
  else:slots=[1,0,1]
  for (op,pos) in reversed(list(zip(word,slots))):
   new={}
   for t,coeff in states.items():
    for pair,v0 in (Q if op=='Q' else H)(t[pos],t[pos+1]).items():
     tt=list(t);tt[pos:pos+2]=pair;tt=tuple(tt)
     new[tt]=new.get(tt,0)+coeff*v0
   states=new
  for t,val in states.items():result[t]=S.expand(result.get(t,0)+sgn*val)
 return {k:S.factor(v) for k,v in result.items()}
for inp,out in [((1,1,2),(1,1,2)),((2,1,1),(2,1,1)),((2,2,1),(2,2,1))]:
 print(inp,out,L(inp,two_Q,two_H).get(out,0))

# triples: 1<2<3
uvars={(a,b):S.symbols(f'u{a}{b}') for a in range(1,4) for b in range(a+1,4)}
tvars={(a,b):S.symbols(f't{a}{b}') for a in range(1,4) for b in range(a+1,4)}
x=S.symbols('x')
def three_Q(a,b):
 if a==b:return {(a,b):l}
 return {(b,a):S.Integer(1),**({(a,b):delta} if a>b else {})}
def three_H(a,b):
 if a==b:return {(a,b):x}
 lo,hi=sorted([a,b]);uv=uvars[lo,hi];tv=tvars[lo,hi]
 return {(b,a):uv} if a<b else {(b,a):tv-uv,(a,b):(q*q+1)*x-q*tv}
for inp,out in [((3,1,2),(3,1,2)),((2,3,1),(2,3,1))]:
 print(inp,out,L(inp,three_Q,three_H).get(out,0))
