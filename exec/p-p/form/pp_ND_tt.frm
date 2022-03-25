off statistics;

s mN,mD,m,t;
in l1,l2,l,mu,nu;
v k1,k2,k3,k4,k;
f pro12,pro32;

l MSQRt = g5_(l1)*(g_(l1,k3)-g_(l1,k1)) * pro12(l1,k1,mN) * g5_(l1)*(g_(l1,k3)-g_(l1,k1)) * pro12(l1,k3,mN) *
          pro12(l2,k2,mN) * pro32(l2,k4,mD,mu,nu) * (k4(mu)-k2(mu)) * (k4(nu)-k2(nu));

id pro12(l?,k?,m?) = g_(l,k) + m*gi_(l);
*id pro32(l?,k?,m?,mu?,nu?) = (g_(l,k) + m*gi_(l)) * (d_(mu,nu)*gi_(l) - 1/3*g_(l,mu,nu) 
*       - 2/(3*m^2)*k(mu)*k(nu)*gi_(l) + 1/(3*m)*(k(mu)*g_(l,nu)-k(nu)*g_(l,mu)) );
id pro32(l?,k?,m?,mu?,nu?) = (g_(l,k) + m*gi_(l)) * (d_(mu,nu)*gi_(l) - 1/3*g_(l,mu,nu) 
       - 2/(3*m^2)*k(mu)*k(nu)*gi_(l) );
.sort

trace4,l1;
trace4,l2;
.sort

id k1.k1 = mN^2;
id k2.k2 = mN^2;
id k3.k3 = mN^2;
id k4.k4 = mD^2;
.sort

id k1.k3 = mN^2-t/2;
id k2.k4 = (mN^2+mD^2-t)/2;
.sort

l MSQRthand = 8*mN^2/(3*mD^2) * t * (t^2 - 2*(mD^2+mN^2)*t + (mD^2-mN^2)^2) *
      ((mD+mN)^2 - t);

l diff = MSQRt - MSQRthand;

l MSQRtDM = -8*mN^2 * 1/(3*mD^2) * t * (t - (mD-mN)^2) * ((mD+mN)^2 - t)^2;

l diff2 = MSQRt - MSQRtDM;

br t;

print;
.end