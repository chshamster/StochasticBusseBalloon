function [X,Y] = circ_cov_sampleB(d,Ntilde)
xi=randn(Ntilde,2)*[1; sqrt(-1)];
Z=fft((d.^0.5).*xi)/sqrt(Ntilde);
X=real(Z); Y=imag(Z);
end 
