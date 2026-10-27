"""Exact, simultaneous D4 check of every tangent in the correctly defined family.

The uniform proof of integrability and independence is in the note; this script
specifically audits the earlier faulty q-derivative computation.
"""
from itertools import combinations,product
from pathlib import Path
import json,time
from check_uniform_identities import build_q,linearized_rows,accum,Laurent


def main():
    start=time.perf_counter();n=4;N=2*n;Q=build_q(n)
    H={src:{} for src in Q};labels=set()
    def mu(i):return tuple([int(i==a)-int(N+1-i==a) for a in range(1,n+1)]+[1])
    def wedge(pair,u,v):
        a,b=map(mu,pair);return a[u]*b[v]-a[v]*b[u]
    for src,col in Q.items():
        for dst,z in col.items():
            f={('scale',):z}
            d=Laurent.make({e:c*e[0] for e,c in z.terms},1)
            if d:f['log_q',]=d
            for u,v in combinations(range(n+1),2):
                c=wedge(dst,u,v)-wedge(src,u,v)
                if c:f['omega',u,v]=c*z
            for a in range(1,n):
                c=dst.count(a)-src.count(a)
                if c:f['local',a]=c*z
            H[src][dst]=f;labels.update(f)
    assert len(labels)==(n+1)*(n+2)//2
    count=0
    for src in product(range(1,N+1),repeat=3):
        rows=linearized_rows(src,Q,H)
        bad={dst:f for dst,f in rows.items() if f}
        assert not bad,(src,bad)
        count+=1
    result={'n':n,'method':'exact simultaneous Laurent expansion',
            'family_directions':len(labels),'tensor_inputs':count,
            'nonzero_linearized_coefficients':0,'elapsed_seconds':round(time.perf_counter()-start,3),
            'status':'PASS','scope':'diagnostic; all-rank integrability proved analytically'}
    (Path(__file__).resolve().parents[1]/'results'/'family_check.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(result,flush=True)

if __name__=='__main__':main()
