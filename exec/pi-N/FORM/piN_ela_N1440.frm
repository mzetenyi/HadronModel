vector pNi,pNo,ppii,ppio,ps,p;
function vertexRNpi,pro1h;
index l;
symbol mN,mR,mpi,m;

local MSQR = pro1h(l,pNo,mN) * vertexRNpi(l,-ppio) * pro1h(l,ps,mR) * vertexRNpi(l,ppii) * 
             pro1h(l,pNi,mN) * vertexRNpi(l,ppii) * pro1h(l,ps,mR) * vertexRNpi(l,-ppio);

id pro1h(l?,p?,m?) = g_(l,p) + gi_(l)*m;
id vertexRNpi(l?,p?) = g5_(l) * g_(l,p);
trace4,l;
.sort

id ps = pNi + ppii;
id ppio = pNi + ppii - pNo;

id pNi.pNi = mN^2;
id pNo.pNo = mN^2;
id ppii.ppii = mpi^2;
id ppio.ppio = mpi^2;

bracket pNi.pNo;
print +s;
.sort

symbol EN,Epi,pabs,srt;
id pNi.ppii = EN*Epi + pabs^2;
id EN = (srt^2 + mN^2 - mpi^2)/(2*srt);
id Epi = (srt^2 - mN^2 + mpi^2)/(2*srt);
id pabs^2 = (srt^4 + mN^4 + mpi^4 - 4*srt^2*(mN^2 + mpi^2) - 4*mN^2*mpi^2)/(4*srt^2);
id srt^2 = mR^2;
id srt^-2 = mR^-2;


bracket pNi.pNo;
print +s;
.end