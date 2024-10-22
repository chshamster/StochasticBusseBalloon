function [ f ] = MostUnst(a,k,m,d,L)
%implementation of Appendix A.1
k=2*pi*k/(2*L);
us=a/2-sqrt(a^2-4*m^2)/2;
vs=(a-us)/m;
A=[-1-vs^2,-2*us*vs;vs^2,-m+2*us*vs];
Gamma=d*A(2,2)+1*A(1,1);
lambda=(Gamma-2*d*k^2)/(d+1);
f=lambda^2+((d+1)*k^2-trace(A))*lambda+d*k^4-Gamma*k^2+det(A);
end
