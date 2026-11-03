const N=8,NAMES=['Personal','Infraestructura','Marketing digital','I+D','Capacitación','Seguridad informática','Soporte y operaciones','Legal y cumplimiento'];
const r6=v=>+v.toFixed(6);
const band=(d,c1,c2)=>d.map((di,i)=>d.map((_,j)=>i===j?di:Math.abs(i-j)===1?c1:Math.abs(i-j)===2?c2:0));
const DS=[30,10,3,1,20,5,2,.5],IA=band([10,8,9,7,8,9,7,8],1,.5),IB=[120,90,100,70,85,95,60,80];
const SC={
ideal:{d:'Presupuesto estable: costos similares y acoplamiento débil entre categorías vecinas. A bien condicionada.',A:IA,b:IB,pc:'jacobi'},
estres:{d:'Inflación de costos: la rigidez de cada categoría sube 10 % y las necesidades b aumentan entre 10 % y 30 %. El equilibrio se desplaza.',A:IA.map((r,i)=>r.map((v,j)=>i===j?r6(v*1.1):v)),b:IB.map((v,i)=>r6(v*[1.15,1.12,1.2,1.3,1.1,1.18,1.25,1.15][i])),pc:'jacobi'},
mal:{d:'Las categorías se miden en escalas muy distintas: A = DCD con D = diag(30, 10, 3, 1, 20, 5, 2, 0,5). Mal escalamiento diagonal: se usa Jacobi.',A:band(DS.map(()=>1),.2,.1).map((r,i)=>r.map((v,j)=>r6(v*DS[i]*DS[j]))),b:[30,15,10,8,12,9,7,5].map((v,i)=>r6(v*DS[i])),pc:'jacobi'}};
let cur='ideal';
const $=id=>document.getElementById(id);
const dot=(a,b)=>a.reduce((s,v,i)=>s+v*b[i],0),nrm=a=>Math.sqrt(dot(a,a));
const mv=(A,x)=>A.map(r=>dot(r,x));
const tr=A=>A[0].map((_,j)=>A.map(r=>r[j]));
function chol(A){const n=A.length,L=A.map(()=>Array(n).fill(0));for(let i=0;i<n;i++)for(let j=0;j<=i;j++){let s=A[i][j];for(let k=0;k<j;k++)s-=L[i][k]*L[j][k];if(i===j){if(!(s>1e-12))return null;L[i][i]=Math.sqrt(s)}else L[i][j]=s/L[j][j]}return L}
function fsub(L,b){const y=[];for(let i=0;i<b.length;i++){let s=b[i];for(let k=0;k<i;k++)s-=L[i][k]*y[k];y[i]=s/L[i][i]}return y}
function bsub(L,y){const n=y.length,x=[];for(let i=n-1;i>=0;i--){let s=y[i];for(let k=i+1;k<n;k++)s-=L[k][i]*x[k];x[i]=s/L[i][i]}return x}
const csolve=(L,b)=>bsub(L,fsub(L,b));
function det(M){const a=M.map(r=>r.slice()),n=a.length;let d=1;for(let c=0;c<n;c++){let p=c;for(let r=c+1;r<n;r++)if(Math.abs(a[r][c])>Math.abs(a[p][c]))p=r;if(Math.abs(a[p][c])<1e-14)return 0;if(p!==c){[a[p],a[c]]=[a[c],a[p]];d=-d}d*=a[c][c];for(let r=c+1;r<n;r++){const f=a[r][c]/a[c][c];for(let k=c;k<n;k++)a[r][k]-=f*a[c][k]}}return d}
function eig(A){const n=A.length,a=A.map(r=>r.slice());for(let s=0;s<100;s++){let off=0;for(let i=0;i<n;i++)for(let j=i+1;j<n;j++)off+=a[i][j]**2;if(off<1e-24)break;for(let p=0;p<n-1;p++)for(let q=p+1;q<n;q++){if(Math.abs(a[p][q])<1e-300)continue;const th=(a[q][q]-a[p][p])/(2*a[p][q]),t=Math.sign(th||1)/(Math.abs(th)+Math.sqrt(th*th+1)),c=1/Math.sqrt(t*t+1),sn=t*c;for(let k=0;k<n;k++){const u=a[k][p],v=a[k][q];a[k][p]=c*u-sn*v;a[k][q]=sn*u+c*v}for(let k=0;k<n;k++){const u=a[p][k],v=a[q][k];a[p][k]=c*u-sn*v;a[q][k]=sn*u+c*v}}}return a.map((r,i)=>r[i])}
const kappa=A=>{const e=eig(A);return Math.max(...e)/Math.min(...e)};
function validate(A,b){const n=A.length;if(A.some(r=>r.some(v=>!Number.isFinite(v)))||b.some(v=>!Number.isFinite(v)))return{ok:false,msg:'Hay celdas vacías o no numéricas en A o b.'};
for(let i=0;i<n;i++)for(let j=i+1;j<n;j++)if(Math.abs(A[i][j]-A[j][i])>1e-9)return{ok:false,msg:`A no es simétrica: a${i+1}${j+1} = ${A[i][j]} y a${j+1}${i+1} = ${A[j][i]}.`};
const minors=A.map((_,k)=>det(A.slice(0,k+1).map(r=>r.slice(0,k+1))));
const L=chol(A);if(!L||minors.some(m=>m<=0))return{ok:false,msg:'A es simétrica pero no definida positiva: algún menor principal es ≤ 0. Corrige los coeficientes.',minors};
return{ok:true,L,minors,eigs:eig(A)}}
function precM(A,t){const n=A.length;if(t==='none')return A.map((_,i)=>A.map((_,j)=>+(i===j)));if(t==='jacobi')return A.map((_,i)=>A.map((_,j)=>i===j?A[i][i]:0));
return A.map((_,i)=>A.map((_,j)=>{let s=0;for(let k=0;k<=Math.min(i,j);k++)s+=A[i][k]*A[k][j]/A[k][k];return s}))}
function pcg(A,b,x0,t,tol,nmax){const M=precM(A,t),LM=chol(M),x=x0.slice(),nb=nrm(b)||1;let r=b.map((v,i)=>v-mv(A,x)[i]);let z=csolve(LM,r),p=z.slice(),rz=dot(r,z);const hist=[nrm(r)/nb];let k=0;
while(hist[k]>=tol&&k<nmax){const Ap=mv(A,p),al=rz/dot(p,Ap);for(let i=0;i<x.length;i++){x[i]+=al*p[i];r[i]-=al*Ap[i]}z=csolve(LM,r);const rzn=dot(r,z),be=rzn/rz;rz=rzn;p=z.map((v,i)=>v+be*p[i]);k++;hist.push(nrm(r)/nb)}
return{x,iters:k,hist,conv:hist[k]<tol,res:hist[k]}}
function kPre(A,t){const L=chol(precM(A,t)),Y=tr(A.map(c=>fsub(L,c))),B=tr(Y.map(r=>fsub(L,r)));return kappa(B)}
function analyze(A,b,x0,t,tol,nmax){const v=validate(A,b);if(!v.ok)return{v};const xr=csolve(v.L,b),g=pcg(A,b,x0,t,tol,nmax),c=pcg(A,b,x0,'none',tol,nmax);
return{v,xr,g,c,kA:kappa(A),kM:kPre(A,t)}}
const f=v=>!Number.isFinite(v)?'—':(v!==0&&(Math.abs(v)>=1e5||Math.abs(v)<1e-4))?v.toExponential(3):(+v.toFixed(5)).toString();
const vec=x=>'['+x.map(f).join(', ')+']';
const sub=['₁','₂','₃','₄','₅','₆','₇','₈'];
function build(){const g=$('grid');g.innerHTML='';for(let i=0;i<N;i++){for(let j=0;j<N;j++)g.insertAdjacentHTML('beforeend',`<input type="number" step="any" id="a${i}${j}" aria-label="a${i+1}${j+1}">`);g.insertAdjacentHTML('beforeend',`<span>|</span><input type="number" step="any" id="b${i}" aria-label="b${i+1}">`)}
g.querySelectorAll('input').forEach(e=>e.addEventListener('input',render));['x0','tol','nmax','pc'].forEach(i=>$(i).addEventListener('input',render))}
function load(k){cur=k;const s=SC[k];s.A.forEach((r,i)=>{r.forEach((v,j)=>$('a'+i+j).value=v);$('b'+i).value=s.b[i]});$('pc').value=s.pc;$('scdesc').textContent=s.d;document.querySelectorAll('[data-sc]').forEach(b=>b.classList.toggle('on',b.dataset.sc===k));render()}
function read(){const A=[],b=[];for(let i=0;i<N;i++){A.push([]);for(let j=0;j<N;j++)A[i].push(parseFloat($('a'+i+j).value));b.push(parseFloat($('b'+i).value))}
let x0=$('x0').value.split(',').map(s=>parseFloat(s));if(x0.length!==N||x0.some(v=>!Number.isFinite(v)))x0=Array(N).fill(0);
let tol=parseFloat($('tol').value);if(!(tol>0))tol=1e-6;let nmax=parseInt($('nmax').value);if(!(nmax>0))nmax=1000;return{A,b,x0,tol,nmax,pc:$('pc').value}}
function table(h,rows){return '<table><tr>'+h.map((x,i)=>`<th${i?'':' class="l"'}>${x}</th>`).join('')+'</tr>'+rows.map(r=>'<tr>'+r.map(c=>`<td>${c}</td>`).join('')+'</tr>').join('')+'</table>'}
function chart(g,c,tol){const W=600,H=220,m=34,mx=Math.max(g.hist.length,c.hist.length)-1||1,all=[...g.hist,...c.hist].map(v=>Math.log10(Math.max(v,1e-16))),lo=Math.floor(Math.min(...all,Math.log10(tol))),hi=Math.ceil(Math.max(...all,0));
const X=i=>m+(W-m-8)*i/mx,Y=v=>8+(H-m-8)*(hi-Math.log10(Math.max(v,1e-16)))/((hi-lo)||1);
const pl=(h,col)=>`<polyline fill="none" stroke="${col}" stroke-width="2" points="${h.map((v,i)=>X(i)+','+Y(v)).join(' ')}"/>`;
let ax='';for(let e=lo;e<=hi;e+=Math.max(1,Math.ceil((hi-lo)/6)))ax+=`<line x1="${m}" x2="${W-8}" y1="${Y(10**e)}" y2="${Y(10**e)}" stroke="var(--ln)"/><text x="2" y="${Y(10**e)+4}">1e${e}</text>`;
return `<svg viewBox="0 0 ${W} ${H}" width="100%" role="img" aria-label="Convergencia del residuo relativo">${ax}<line x1="${m}" x2="${W-8}" y1="${Y(tol)}" y2="${Y(tol)}" stroke="var(--er)" stroke-dasharray="4 3"/>${pl(c.hist,'var(--wa)')}${pl(g.hist,'var(--ac)')}<text x="${m}" y="${H-10}">iteración 0 → ${mx}</text><text x="${W-210}" y="14" style="fill:var(--wa)">GC simple</text><text x="${W-140}" y="14" style="fill:var(--ac)">GCP</text><text x="${W-90}" y="14" style="fill:var(--er)">tolerancia</text></svg>`}
function render(){setTimeout(qa,0);const p=read(),R=analyze(p.A,p.b,p.x0,p.pc,p.tol,p.nmax);
$('sys').textContent=p.A.map((r,i)=>r.map((v,j)=>(j?(v<0?' − ':' + '):'')+f(j?Math.abs(v):v)+'·x'+sub[j]).join('')+'  =  '+f(p.b[i])).join('\n');
if(!R.v.ok){$('alert').innerHTML=`<div class="msg" role="alert"><strong>Matriz no válida.</strong> ${R.v.msg} No se calcularon resultados.</div>`;$('verif').innerHTML=`<div class="msg">${R.v.msg}</div>`;['ref','cards','gcpt','chart','err','errx'].forEach(i=>$(i).innerHTML='');cmp(p);return}
$('alert').innerHTML='<div class="msg good">Matriz válida: simétrica y definida positiva.</div>';
const v=R.v;$('verif').innerHTML=`<p>Simetría: aᵢⱼ = aⱼᵢ para todo i, j. <span class="ok">Cumple.</span></p>`+table(['Criterio','Valores','Resultado'],[['Menores principales (Sylvester)',v.minors.map(f).join(', '),'<span class="ok">todos &gt; 0</span>'],['Autovalores',v.eigs.map(f).join(', '),'<span class="ok">todos &gt; 0</span>']]).replace(/<td>(?=[^<]*<\/td><td>)/,'<td class="l">');
$('ref').innerHTML=table(['Variable','Categoría','x de referencia (miles de Bs)'],R.xr.map((x,i)=>['x'+sub[i],NAMES[i],f(x)]));
const pn={jacobi:'Jacobi',ssor:'SSOR',none:'ninguno'}[p.pc];
$('cards').innerHTML=[['GC simple',R.c],['GCP ('+pn+')',R.g]].map(([n,r])=>`<div><div class="ledger">${n}</div><div class="big">${r.iters} iter.</div><div>${r.conv?'<span class="ok">Convergió</span>':'<span class="bad">No convergió</span>'} · residuo ${f(r.res)}</div></div>`).join('')+`<div><div class="ledger">Condicionamiento</div><div>κ(A) = <strong>${f(R.kA)}</strong></div><div>κ(M⁻¹A) = <strong>${f(R.kM)}</strong></div></div>`;
$('gcpt').innerHTML=table(['Método','Vector aproximado','Iteraciones','Residuo relativo final','¿Convergió?'],[['GC simple',vec(R.c.x),R.c.iters,f(R.c.res),R.c.conv?'Sí':'No'],['GCP',vec(R.g.x),R.g.iters,f(R.g.res),R.g.conv?'Sí':'No']]);
$('chart').innerHTML=chart(R.g,R.c,p.tol);
$('err').innerHTML=table(['Variable','x referencia','x GCP','Error absoluto','Error relativo','Error porcentual (%)'],R.xr.map((x,i)=>{const ea=Math.abs(x-R.g.x[i]);return['x'+sub[i],f(x),f(R.g.x[i]),f(ea),x===0?'no definido':f(ea/Math.abs(x)),x===0?'no definido':f(ea/Math.abs(x)*100)]}));
const rr=nrm(p.b.map((v,i)=>v-mv(p.A,R.g.x)[i]))/nrm(p.b);$('errx').innerHTML=`<p>Residuo relativo global ‖b − Ax_GCP‖₂ / ‖b‖₂ = <strong>${f(rr)}</strong></p>`;cmp(p)}
const INT={ideal:'Equilibrio estable: la asignación refleja las necesidades y converge con pocas iteraciones.',estres:'La inflación desplaza el equilibrio hacia montos mayores; el condicionamiento casi no cambia.',mal:'El mal escalamiento eleva κ(A); el precondicionador recupera una convergencia rápida y estable.'};
function cmp(p){const names={ideal:'Ideal',estres:'Bajo estrés',mal:'Mal condicionado'},R={};for(const k in SC)R[k]=analyze(SC[k].A,SC[k].b,p.x0,p.pc,p.tol,p.nmax);
const row=(l,fn)=>[l,...Object.keys(SC).map(k=>R[k].v.ok?fn(R[k]):'—')];
$('cmp').innerHTML=table(['Indicador',...Object.values(names)],[
row('Solución de referencia',r=>vec(r.xr)),row('Solución GCP',r=>vec(r.g.x)),
row('Error absoluto (Ea)',r=>vec(r.xr.map((x,i)=>Math.abs(x-r.g.x[i])))),
row('Error porcentual (Ep %)',r=>vec(r.xr.map((x,i)=>x===0?NaN:Math.abs(x-r.g.x[i])/Math.abs(x)*100))),
row('Iteraciones GC',r=>r.c.iters),row('Iteraciones GCP',r=>r.g.iters),row('κ(A)',r=>f(r.kA)),row('κ(M⁻¹A)',r=>f(r.kM)),row('Residuo relativo final GCP',r=>f(r.g.res)),row('¿Convergió?',r=>r.g.conv?'Sí':'No'),
['Interpretación financiera',...Object.keys(SC).map(k=>`<span style="white-space:normal;display:block;min-width:200px;text-align:left">${INT[k]}</span>`)]])}
function qa(){const p=read(),S=k=>SC[k],T=[],P={x0:Array(N).fill(0),tol:1e-6,nmax:1000};
const run=k=>analyze(S(k).A,S(k).b,P.x0,'jacobi',P.tol,P.nmax);
let o=['ideal','estres','mal'].map(k=>{const r=run(k);return k+': '+(r.v.ok&&r.g.conv?'OK ('+r.g.iters+' it.)':'falla')});
T.push(['Carga de escenarios','Ideal, Bajo estrés y Mal condicionado','Los 3 son válidos y GCP converge',o.join('; '),o.every(s=>s.includes('OK'))]);
const A1=S('ideal').A.map(r=>r.slice());A1[0][1]=5;const v1=validate(A1,S('ideal').b);
const A2=Array.from({length:N},(_,i)=>Array.from({length:N},(_,j)=>+(i===j)));A2[0][1]=A2[1][0]=3;const v2=validate(A2,Array(N).fill(1));
T.push(['Matriz no válida','Caso 1: a₁₂=5, a₂₁=1. Caso 2: simétrica con menor principal negativo','Ambos casos se rechazan con mensaje, sin colapsar',`Caso 1: ${v1.ok?'aceptada':v1.msg} | Caso 2: ${v2.ok?'aceptada':v2.msg}`,!v1.ok&&!v2.ok]);
const A3=S('ideal').A.map(r=>r.slice());A3[0][0]=20;const r0=run('ideal'),r3=analyze(A3,S('ideal').b,P.x0,'jacobi',P.tol,P.nmax);
const diff=r3.v.ok&&Math.abs(r3.xr[0]-r0.xr[0])>1e-3,close=r3.v.ok&&nrm(r3.xr.map((x,i)=>x-r3.g.x[i]))<1e-4;
T.push(['Edición dinámica','a₁₁: 10 → 20 en el escenario ideal','x₁ cambia y GCP coincide con la referencia',`x₁: ${f(r0.xr[0])} → ${f(r3.xr[0])}; ‖xref − xGCP‖ ${f(nrm(r3.xr.map((x,i)=>x-r3.g.x[i])))}`,diff&&close]);
const rm=run('mal');
T.push(['Comparativa de precondicionamiento','Escenario mal condicionado, Jacobi','Iteraciones GCP < iteraciones GC (ganancia real)',`GCP ${rm.g.iters} < GC ${rm.c.iters}; κ(A)=${f(rm.kA)} → κ(M⁻¹A)=${f(rm.kM)}`,rm.g.iters<rm.c.iters]);
const vc=validate(p.A,p.b);T.push(['Validación de datos actuales','Matriz A y vector b cargados ahora en la interfaz','A simétrica y definida positiva',vc.ok?'Matriz válida':'La matriz actual no es SPD o no es simétrica. '+vc.msg,vc.ok,'REVISAR']);
$('qat').innerHTML=table(['Prueba','Entrada','Resultado esperado','Resultado obtenido','Conclusión'],T.map(t=>[t[0],t[1],t[2],t[3],t[4]?'<span class="ok">Aprobada</span>':'<span class="bad">'+(t[5]||'Fallida')+'</span>'])).replace(/<td>/g,'<td class="w">')}
build();document.querySelectorAll('[data-sc]').forEach(b=>b.addEventListener('click',()=>load(b.dataset.sc)));$('runqa').addEventListener('click',qa);load('ideal');qa();
