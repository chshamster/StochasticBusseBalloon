function [c] = Cov(x,ell)
c=1/(2*ell)*exp(-pi*x.^2/(4*ell^2));
end

