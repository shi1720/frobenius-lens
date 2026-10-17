"""Standalone verifier primitives: no SymPy or certificate code imports."""
def evolve(rule, n, t, seed, r, p):
    u = list(seed) + [0]*(n-len(seed))
    active = [(j-r, c) for j,c in enumerate(rule) if c % p]
    for _ in range(t):
        u = [sum(c*u[(i-shift)%n] for shift,c in active)%p for i in range(n)]
    return tuple(u)

def rank(A, p):
    A = [list(row) for row in A]
    if not A:
        return 0
    m,n = len(A),len(A[0]); k=0
    for j in range(n):
        pivot=next((i for i in range(k,m) if A[i][j]%p),None)
        if pivot is None:
            continue
        A[k],A[pivot] = A[pivot],A[k]
        inv=pow(A[k][j],-1,p)
        A[k]=[(v*inv)%p for v in A[k]]
        for i in range(k+1,m):
            a=A[i][j]
            if a:
                A[i]=[(v-a*w)%p for v,w in zip(A[i],A[k])]
        k+=1
        if k==m:
            break
    return k

def observation_matrix(specs, r, p):
    cols=[]
    for j in range(2*r+1):
        f=[0]*(2*r+1); f[j]=1
        cols.append(sum((evolve(f,n,t,s,r,p) for n,t,s in specs),()))
    return list(zip(*cols))

def digits(i,n,p):
    out=[]
    for _ in range(n):
        out.append(i%p); i//=p
    return tuple(out)
