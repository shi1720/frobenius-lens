"""Execute the frozen protocol from the repository root."""
import sys, json, csv, hashlib, platform, time, random, statistics
from pathlib import Path
from collections import Counter
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'src'), str(ROOT/'tests')]
from frobenius_lens import Experiment, certificate, snapshot, minimum_repair, annihilator
from independent import evolve, rank, observation_matrix, digits
P=json.loads((ROOT/'experiments/protocol.json').read_text())
OUT=ROOT/'results'; OUT.mkdir(exist_ok=True)
start=time.perf_counter(); rng=random.Random(P['suite_rank_checks']['rng_seed'])
counts=Counter(); failures=[]

def check(label, condition, info):
    counts[label]+=1
    if not condition:
        failures.append({'check':label,'info':info})
        raise AssertionError(failures[-1])

with (OUT/'single_checks.csv').open('w') as f:
    w=csv.writer(f); w.writerow(['p','n','t','seed_integer','radius','degree','predicted_rank','direct_rank'])
    for conf in P['single']:
        p=conf['p']
        for n in range(1,conf['n_max']+1):
            for a in range(conf['a_max']+1):
                t=p**a
                for seed_id in range(p**n):
                    s=digits(seed_id,n,p); e=Experiment(n,t,s)
                    for r in range(conf['r_max']+1):
                        c=certificate([e],r,p)
                        k=rank(observation_matrix([(n,t,s)],r,p),p)
                        check('single_rank',k==c['rank'],[p,n,t,seed_id,r])
                        w.writerow([p,n,t,seed_id,r,c['degree'],c['rank'],k])
        print('single field complete',p,dict(counts),flush=True)

suite_rows=[]
for p,N in [(2,400),(3,200)]:
    for trial in range(N):
        r=rng.randrange(9); es=[]
        for _ in range(rng.randrange(1,5)):
            n=rng.randrange(1,25); t=p**rng.randrange(4)
            s=tuple(rng.randrange(p) for _ in range(n))
            es.append(Experiment(n,t,s))
        c=certificate(es,r,p)
        k=rank(observation_matrix([(e.n,e.t,e.seed) for e in es],r,p),p)
        check('suite_rank',k==c['rank'],[p,trial])
        rule=tuple(rng.randrange(p) for _ in range(2*r+1))
        for e in es:
            check('snapshot_agreement',snapshot(rule,e,r,p)==evolve(rule,e.n,e.t,e.seed,r,p),[p,trial])
        for h in c['kernel_basis']:
            h=tuple(h+[0]*(2*r+1-len(h)))
            for e in es:
                check('kernel_witness',all(v==0 for v in evolve(h,e.n,e.t,e.seed,r,p)),[p,trial])
        suite_rows.append({'p':p,'trial':trial,'radius':r,'certificate':c})
(OUT/'suite_checks.json').write_text(json.dumps(suite_rows,separators=(',',':')))
print('suite checks complete',dict(counts),flush=True)

fibers=[]
for p,radii,nmax,times in [(2,[1,2,3],6,[1,2,4,8]),(3,[1],5,[1,3,9])]:
    for r in radii:
        for n in range(1,nmax+1):
            seeds=[(0,)*n,(1,)+(0,)*(n-1),(1,)*n,tuple(i%p for i in range(n)),tuple(rng.randrange(p) for _ in range(n))]
            for t in times:
                for s in seeds:
                    e=Experiment(n,t,s); c=certificate([e],r,p)
                    hist=Counter(evolve(digits(i,2*r+1,p),n,t,s,r,p) for i in range(p**(2*r+1)))
                    check('fiber_enumeration',len(hist)==p**c['rank'] and set(hist.values())=={c['fiber_size']},[p,r,n,t,s])
                    fibers.append({'p':p,'r':r,'n':n,'t':t,'seed':s,'image_size':len(hist),'fiber_size':c['fiber_size']})
(OUT/'fiber_checks.json').write_text(json.dumps(fibers,separators=(',',':')))
print('fiber checks complete',dict(counts),flush=True)

# Analytical impulse formula, independently compared across a wider geometry grid.
impulse=[]
for p in [2,3,5]:
    for n in range(1,129):
        for a in range(6):
            e=Experiment(n,p**a); c=certificate([e],4,p)
            b,u=0,n
            while u%p==0:
                b,u=b+1,u//p
            effective=n//(p**min(a,b))
            check('impulse_formula',c['degree']==effective,[p,n,a])
            impulse.append([p,n,a,effective,c['rank']])
with (OUT/'impulse_grid.csv').open('w') as f:
    w=csv.writer(f); w.writerow(['p','n','a','effective_length','rank_r4']); w.writerows(impulse)

case=P['case_study']; r=case['radius']; p=case['p']; t=case['time']; target=case['target']
base=[Experiment(n,t) for n in case['base_rings']]
strategies={'base':base,'more_dyadic':base+[Experiment(256,t)],'earlier_snapshot':base+[Experiment(128,1)],'odd_repair':base+[Experiment(9,t)],'bounded_repair':base+[Experiment(3,t),Experiment(7,t)]}
case_results={}
for name,es in strategies.items():
    c=certificate(es,r,p)
    ref=[evolve(target,e.n,e.t,e.seed,r,p) for e in es]
    accepted=[]
    for i in range(2**(2*r+1)):
        rule=digits(i,2*r+1,p)
        if all(evolve(rule,e.n,e.t,e.seed,r,p)==y for e,y in zip(es,ref)):
            accepted.append(i)
    check('case_fiber',len(accepted)==c['fiber_size'],name)
    case_results[name]={'certificate':c,'accepted_rule_ids':accepted,'detected_alternatives':512-len(accepted),'missed_alternatives':len(accepted)-1,'stored_cells':sum(e.n for e in es)}
    print('case',name,'fiber',len(accepted),flush=True)
menu=[Experiment(n,t) for n in case['repair_menu']]
repair=minimum_repair(base,r,p,menu)
bounded_menu=[Experiment(n,t) for n in range(1,8)]
bounded=minimum_repair(base,r,p,bounded_menu)
check('repair_full',repair['cost']==9 and repair['indices']==[8],repair)
check('repair_bounded',bounded['cost']==10 and bounded['indices']==[2,6],bounded)
case_results['repair_search']={'full_menu':repair,'bounded_menu':bounded}
(OUT/'case_study.json').write_text(json.dumps(case_results,indent=2))

# Guard against unsupported inputs; no silently generalized certificates.
for bad,p in [(Experiment(8,3),2),(Experiment(0,1),2),(Experiment(8,1,(2,)),2),(Experiment(8,1),4)]:
    try:
        certificate([bad],1,p)
    except ValueError:
        check('domain_guard',True,str(bad))
    else:
        check('domain_guard',False,str(bad))
# t=3 is genuinely outside the linear-in-rule setting.
f,g=(1,1,0),(0,1,1); h=tuple((a+b)%2 for a,b in zip(f,g))
e=Experiment(7,1)
yf=evolve(f,7,3,(1,),1,2); yg=evolve(g,7,3,(1,),1,2); yh=evolve(h,7,3,(1,),1,2)
check('non_frobenius_counterexample',yh!=tuple((a+b)%2 for a,b in zip(yf,yg)),[yf,yg,yh])

# Descriptive timing: cached certificates and independent basis simulation, separate tasks.
timings=[]
for n in [32,64,128,256]:
    e=Experiment(n,64); certificate([e],4,2)
    a=[];b=[]
    for _ in range(5):
        z=time.perf_counter(); certificate([e],4,2);a.append(time.perf_counter()-z)
        z=time.perf_counter(); rank(observation_matrix([(n,64,(1,))],4,2),2);b.append(time.perf_counter()-z)
    timings.append({'n':n,'certificate_cached_seconds':a,'direct_basis_seconds':b})
(OUT/'timings.json').write_text(json.dumps(timings,indent=2))
summary={'counts':dict(counts),'total_checks':sum(counts.values()),'failures':failures,'wall_seconds':time.perf_counter()-start,'protocol_sha256':hashlib.sha256((ROOT/'experiments/protocol.json').read_bytes()).hexdigest(),'python':platform.python_version(),'platform':platform.platform(),'notes':'Exact tests; seeds are deterministic. No population error-rate or LLM performance inference.'}
(OUT/'summary.json').write_text(json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2))
