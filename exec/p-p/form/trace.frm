off statistics;

vector k1,k2,k3;
index l,mu;

local T1 = g_(l,k1,k3,k2,mu);

local T = g_(l,k1,k3,k2,mu) - g_(l,k2,k3,k1,mu);

trace4,l;

print;

.end

