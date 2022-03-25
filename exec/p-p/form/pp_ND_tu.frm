off statistics;

s mN,mD,m;
s s,t,u;
in l;
f P,pro12;
v k1,k2,k3,k4,pt,pu,k,p1,p2;

l MSQRut = 
+ g5_(l)*pro12(l,k1,mN)*pro12(l,k4,mD)*P(l,k4,pu,pt)*pro12(l,k2,mN)*g5_(l)*pro12(l,k3,mN)
+ g5_(l)*pro12(l,k2,mN)*pro12(l,k4,mD)*P(l,k4,pt,pu)*pro12(l,k1,mN)*g5_(l)*pro12(l,k3,mN);

id pro12(l?,k?,m?) = g_(l,k) + gi_(l)*m;
id P(l?,k?,p1?,p2?) = - gi_(l)*p1.p2 + 1/3*g_(l,p1,p2) + 2/(3*mD^2)*gi_(l)*k.p1*k.p2 - 1/(3*mD)*(k.p1*g_(l,p2)-k.p2*g_(l,p1));

.sort

trace4,l;
.sort

id pt = k1 - k3;
id pu = k2 - k3;
id k4 = k1 + k2 - k3;

id k1.k1 = mN^2;
id k2.k2 = mN^2;
id k3.k3 = mN^2;

id k1.k2 = s/2 - mN^2;
id k1.k3 = mN^2 - t/2;
id k2.k3 = mN^2 - u/2;
id s = 3*mN^2 + mD^2 - t - u;
.sort

l MSQRutDM1 = 1/(2*mD^2) * (
    (t*u + (mD^2-mN^2)*(t+u) - mD^4 + mN^4) * (t*u + mN*(mN+mD)*(mD^2-mN^2))
   + 1/3 * (t*u - (mD+mN)^2*(t+u) + (mN+mD)^4) * (t*u - mN*(mD-mN)*(mD^2-mN^2))
);

l MSQRutDM2 = 1/(2*mD^2) * (
    (t*u + 2*(mD^2-mN^2)*(t+u) - mD^4 + mN^4) * (t*u + mN*(mN+mD)*(mD^2-mN^2))
   + 1/3 * (t*u + 2*(mD+mN)^2*(t+u) + (mN+mD)^4) * (t*u - mN*(mD-mN)*(mD^2-mN^2))
);

l MSQRutDM3 = 1/(2*mD^2) * (
    (t*u + 2*(mD^2-mN^2)*(t+u) - (mD^2-mN^2)*(2*mD^2-mD*mN+3*mN^2)) * (t*u + mN*(mN+mD)*(mD^2-mN^2))
   + 1/3 * (t*u + 2*(mD+mN)^2*(t+u) - (mD+mN)^2*(2*mD^2 + mD*mN + 5*mN^2)) * (t*u - mN*(mD-mN)*(mD^2-mN^2))
);

l diff1 = MSQRut - MSQRutDM1;
l diff2 = MSQRut - MSQRutDM2;
l diff3 = MSQRut - MSQRutDM3;

br t,u;

print;

.end