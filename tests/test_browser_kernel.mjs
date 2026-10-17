import fs from 'node:fs';
import assert from 'node:assert/strict';
import {annihilator,degree,certificate,snapshot,asCoefficients,step} from '../app/lens.mjs';
let count=0;
for(let line of fs.readFileSync(new URL('../results/single_checks.csv',import.meta.url),'utf8').trim().split('\n').slice(1)){
 let [p,n,t,s,r,D,rank]=line.split(',').map(Number);if(p!==2)continue;
 let e={n,t,seed:BigInt(s)};
 assert.equal(degree(annihilator(e)),D);assert.equal(certificate([e],r).rank,rank);count++;
}
let es=[32,64,128].map(n=>({n,t:64,seed:1n}));
assert.equal(certificate(es,4).fiber,128n);
assert.equal(certificate([...es,{n:9,t:64,seed:1n}],4).fiber,1n);
for(let n=1;n<=12;n++)for(let t of [1,2,4,8])for(let i=0;i<32;i++){
 let f=asCoefficients(BigInt(i),5),e={n,t,seed:1n},u=asCoefficients(1n,n);
 for(let j=0;j<t;j++)u=step(u,f,2);
 assert.deepEqual(asCoefficients(snapshot(f,e,2),n),u);count++;
}
fs.writeFileSync(new URL('../results/browser_kernel_checks.json',import.meta.url),JSON.stringify({checks:count+2,failures:0,description:'Binary JS kernel against Python degree/rank fixtures plus independent repeated time stepping.'},null,2));
console.log(`${count+2} browser kernel checks passed`);
