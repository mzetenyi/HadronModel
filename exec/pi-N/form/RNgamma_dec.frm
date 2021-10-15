v pi,pf,k;
in mu,mup,al,alp;
in l;
s mN,mR;

local MSQR = -(g_(l,pf)+mN*gi_(l)) * i_/2*(g_(l,k,mu)-g_(l,mu,k)) * (g_(l,pi)+mR*gi_(l)) * i_/2*(g_(l,k,mu)-g_(l,mu,k));

trace4,l;
print;
.end