"""Exact certificates. Coefficients are in increasing degree order throughout."""
from dataclasses import dataclass
from functools import lru_cache
from itertools import combinations
import sympy as sp

X = sp.Symbol('x')

def prime(p):
    if not isinstance(p, int) or not sp.isprime(p):
        raise ValueError('p must be prime')

def poly(coeff, p):
    return sp.Poly.from_list(list(reversed(coeff)) or [0], X, modulus=p)

def coeff(q, p):
    return [int(q.nth(i)) % p for i in range(int(q.degree()) + 1)] if q else [0]

def frobenius_exponent(t, p):
    if not isinstance(t, int) or t < 1:
        raise ValueError('time must be a positive power of p (including 1)')
    a, u = 0, t
    while u % p == 0:
        a, u = a + 1, u // p
    if u != 1:
        raise ValueError('certificate supports only times p**a')
    return a

@dataclass(frozen=True)
class Experiment:
    n: int
    t: int
    seed: tuple[int, ...] = (1,)

    def validate(self, p):
        prime(p)
        if not isinstance(self.n, int) or self.n < 1:
            raise ValueError('ring length must be positive')
        frobenius_exponent(self.t, p)
        if len(self.seed) > self.n or any(not isinstance(c, int) or not 0 <= c < p for c in self.seed):
            raise ValueError('seed must have at most n canonical field coefficients')

@lru_cache(maxsize=4096)
def annihilator(e, p):
    """Q_e divides h iff this experiment cannot detect rule difference x^-r h."""
    e.validate(p)
    m = sp.Poly(X**e.n - 1, X, modulus=p)
    s = poly(e.seed, p)
    q = m.exquo(sp.gcd(m, s))
    out = sp.Poly(1, X, modulus=p)
    for f, v in sp.factor_list(q)[1]:
        out *= f ** ((int(v) + e.t - 1) // e.t)
    return out.monic()

def certificate(experiments, radius, p=2):
    prime(p)
    if not isinstance(radius, int) or radius < 0:
        raise ValueError('radius must be a nonnegative integer')
    L = sp.Poly(1, X, modulus=p)
    details = []
    for e in experiments:
        Q = annihilator(e, p)
        L = sp.lcm(L, Q).monic()
        details.append({'n': e.n, 't': e.t, 'seed': list(e.seed), 'Q': coeff(Q,p)})
    d, D = 2 * radius, int(L.degree())
    nullity = max(0, d + 1 - D)
    basis = [coeff(L * sp.Poly(X**j, X, modulus=p), p) for j in range(nullity)]
    return {'p': p, 'radius': radius, 'experiments': details,
            'L': coeff(L,p), 'degree': D, 'rank': d+1-nullity,
            'nullity': nullity, 'fiber_size': p**nullity,
            'identifiable': nullity == 0, 'kernel_basis': basis}

def snapshot(rule, e, radius, p=2):
    """Compute a snapshot by polynomial powering, for comparison with time stepping."""
    e.validate(p)
    if len(rule) != 2*radius+1:
        raise ValueError('rule must have 2r+1 coefficients')
    M = sp.Poly(X**e.n - 1, X, modulus=p)
    # x is a unit modulo M. Use a nonnegative representative of x^-r.
    f = sp.rem(poly(rule,p) * sp.Poly(X**((-radius) % e.n), X, modulus=p), M)
    y, b, t = sp.Poly(1,X,modulus=p), f, e.t
    while t:
        if t & 1:
            y = sp.rem(y*b, M)
        b, t = sp.rem(b*b,M), t//2
    y = sp.rem(y*poly(e.seed,p),M)
    return tuple(int(y.nth(i)) % p for i in range(e.n))

def minimum_repair(experiments, radius, p, candidates, costs=None):
    """Exact finite-menu search; exponential in menu size, no greedy claim."""
    if len(candidates) > 22:
        raise ValueError('exact subset search is capped at 22 candidates')
    costs = [e.n for e in candidates] if costs is None else list(costs)
    if len(costs) != len(candidates) or any(c <= 0 for c in costs):
        raise ValueError('one positive cost is required per candidate')
    best = None
    base = sp.Poly(1, X, modulus=p)
    for e in experiments:
        base = sp.lcm(base, annihilator(e,p)).monic()
    if int(base.degree()) > 2*radius:
        return {'indices': [], 'cost': 0, 'degree': int(base.degree())}
    qs = [annihilator(e,p) for e in candidates]
    for k in range(1,len(candidates)+1):
        for ids in combinations(range(len(candidates)),k):
            cost = sum(costs[i] for i in ids)
            if best is not None and cost >= best['cost']:
                continue
            q = base
            for i in ids:
                q = sp.lcm(q,qs[i]).monic()
            if int(q.degree()) > 2*radius:
                best = {'indices': list(ids), 'cost': cost, 'degree': int(q.degree())}
    return best

def recover(experiments, observations, radius, p=2):
    """Return one compatible rule and a complete basis of rule ambiguity.

    Reject inconsistent observations. This validates data against the declared
    scalar linear class; it does not certify an arbitrary simulator's source.
    """
    es=list(experiments); ys=list(observations)
    c=certificate(es,radius,p)
    if len(es)!=len(ys) or any(len(y)!=e.n for e,y in zip(es,ys)):
        raise ValueError('one full-ring observation is required per experiment')
    if any(not isinstance(v,int) or not 0<=v<p for y in ys for v in y):
        raise ValueError('observations must use canonical field coefficients')
    k=2*radius+1
    cols=[]
    for j in range(k):
        f=[0]*k;f[j]=1
        cols.append(sum((snapshot(f,e,radius,p) for e in es),()))
    A=[list(row)+[b] for row,b in zip(zip(*cols),sum((tuple(y) for y in ys),()))]
    pivots=[]; row=0
    for j in range(k):
        pivot=next((i for i in range(row,len(A)) if A[i][j]%p),None)
        if pivot is None:
            continue
        A[row],A[pivot]=A[pivot],A[row]
        z=pow(A[row][j],-1,p);A[row]=[(z*v)%p for v in A[row]]
        for i in range(len(A)):
            if i!=row and A[i][j]:
                z=A[i][j];A[i]=[(v-z*w)%p for v,w in zip(A[i],A[row])]
        pivots.append(j);row+=1
    if any(not any(a[:k]) and a[k] for a in A):
        raise ValueError('observations are inconsistent with the declared rule class')
    particular=[0]*k
    for i,j in enumerate(pivots):
        particular[j]=A[i][-1]
    if row!=c['rank']:
        raise ArithmeticError('symbolic certificate and observation rank disagree')
    return {'particular_rule':particular,'certificate':c}
