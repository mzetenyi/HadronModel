f vertexRNpi, vertexRNpiAdj;
f vertexRNgamma, vertexRNgammaAdj;
f proN, propR;
f T,Tadj;
s mR,mN,mpi;
s srt;
v pi,k,p,pf,q;
in l,mu;

local MSQR = proN(l,pf)*T(l,pi,k,p,pf,q,mu)*proN(l,pi)*Tadj(l,pi,k,p,pf,q,mu);

id T(l?,pi?,k?,p?,pf?,q?,mu?) = vertexRNpi(l,p,pf,q) * propR(l,p) * vertexRNgamma(l,p,pi,k,mu);
id Tadj(l?,pi?,k?,p?,pf?,q?,mu?) = vertexRNgammaAdj(l,p,pi,k,mu) * propR(l,p) * vertexRNpiAdj(l,p,pf,q);

id vertexRNpi(l?,p?,pf?,q?) = g5_(l) * g_(l,q);
id vertexRNpiAdj(l?,p?,pf?,q?) = g5_(l) * g_(l,q);
id vertexRNgamma(l?,p?,pi?,k?,mu?) = g_(l,k,mu)-g_(l,mu,k);
id vertexRNgammaAdj(l?,p?,pi?,k?,mu?) = -(g_(l,k,mu)-g_(l,mu,k));

id proN(l?,p?) = g_(l,p) + mN*gi_(l);
id propR(l?,p?) = g_(l,p) + mR*gi_(l);

trace4,l;
.sort

id k.k=0;
id pi.pi=mN^2;
id pf.pf=mN^2;
id q.q=mpi^2;
id p.p=srt^2;

print +s;
.end