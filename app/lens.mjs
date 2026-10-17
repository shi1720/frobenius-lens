// Binary, factorization-free implementation independent of the Python backend.
export const degree = a => a===0n ? -1 : a.toString(2).length-1;
export function mod(a,b){if(b===0n)throw Error('zero divisor');let db=degree(b);while(degree(a)>=db)a^=b<<BigInt(degree(a)-db);return a;}
export function mul(a,b){let z=0n;while(b){if(b&1n)z^=a;a<<=1n;b>>=1n;}return z;}
export function gcd(a,b){while(b){[a,b]=[b,mod(a,b)];}return a;}
export function quotient(a,b){let q=0n,db=degree(b);if(db<0)throw Error('zero divisor');while(degree(a)>=db){let s=BigInt(degree(a)-db);q^=1n<<s;a^=b<<s;}if(a)throw Error('non-exact division');return q;}
export const lcm=(a,b)=>mul(quotient(a,gcd(a,b)),b);
export function rotate(a,k,n){k=((k%n)+n)%n;const mask=(1n<<BigInt(n))-1n;return ((a<<BigInt(k))|(a>>BigInt(n-k)))&mask;}
export function seedBits(n,mode){if(mode==='zero')return 0n;if(mode==='constant')return (1n<<BigInt(n))-1n;if(mode==='impulse')return 1n;let s=0n,z=(1720406^n)>>>0;for(let i=0;i<n;i++){z^=z<<13;z^=z>>>17;z^=z<<5;if(z&1)s|=1n<<BigInt(i);}return s;}
export function validate(e){if(!Number.isInteger(e.n)||e.n<1||e.n>512)throw Error('Use ring lengths 1–512.');if(!Number.isInteger(e.t)||e.t<1||e.t>1048576||(e.t&(e.t-1)))throw Error('Time must be a power of two, at most 2^20.');if(e.seed<0n||degree(e.seed)>=e.n)throw Error('Seed is outside the ring.');}
export function annihilator(e){validate(e);let basis=new Map();for(let j=0;j<=e.n;j++){let v=rotate(e.seed,(j*e.t)%e.n,e.n),c=1n<<BigInt(j);while(v){let k=degree(v),b=basis.get(k);if(!b){basis.set(k,[v,c]);break;}v^=b[0];c^=b[1];}if(!v)return c;}throw Error('No dependency found');}
export function certificate(es,r){if(!Number.isInteger(r)||r<0||r>32)throw Error('Radius must be 0–32.');let L=1n;let rows=es.map(e=>{let Q=annihilator(e);L=lcm(L,Q);return {...e,Q};});let nullity=Math.max(0,2*r+1-degree(L));return {L,degree:degree(L),rank:2*r+1-nullity,nullity,fiber:1n<<BigInt(nullity),rows};}
export function smallestRepair(es,r,max=64){let c=certificate(es,r);if(!c.nullity)return null;for(let n=1;n<=max;n++){let e={n,t:es[0]?.t||64,seed:1n,mode:'impulse'};if(degree(lcm(c.L,annihilator(e)))>2*r)return e;}return null;}
export function step(u,rule,r){let n=u.length;return u.map((_,i)=>rule.reduce((z,c,j)=>z^(c&u[((i-j+r)%n+n)%n]),0));}
export function snapshot(rule,e,r){validate(e);let out=0n;rule.forEach((c,j)=>{if(c&1)out^=rotate(e.seed,((j-r)*e.t)%e.n,e.n);});return out;}
export function asCoefficients(a,n=degree(a)+1){return Array.from({length:n},(_,i)=>Number((a>>BigInt(i))&1n));}
export function polynomial(a){if(a===0n)return '0';return asCoefficients(a).flatMap((v,i)=>v?[i===0?'1':i===1?'x':`x^${i}`]:[]).reverse().join(' + ');}
