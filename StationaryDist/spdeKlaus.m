function [PulseCount]=spdeKlaus(IC,t,dt,L,x,J,h,a,sigma,q,dEE,S1,MAX,wl,Start,End)
u_n=IC; % set initial condition
sqrtdts=sqrt(dt)*sigma;
%compute settings for Lord's algorithms
tilde_q=[q; q(end-1:-1:2)];
Ntilde=length(tilde_q);
d=ifft(tilde_q,'symmetric')*Ntilde;
flag=false;
for k=1:length(t)-1 % time loop
    if flag==false % generate two samples
        [dW11,dW21]=circ_cov_sampleB(d,Ntilde);
        dW1=cat(1,zeros(J+1,1),dW11(1:J+1)); dW2=cat(1,zeros(J+1,1),dW21(1:J+1));
        flag=true;
    else % use second sample from last call
        dW1=dW2; flag=false;
    end
    fu=f(u_n,a,J); % evaluate f
    u_new=dEE\(u_n+dt*fu+sqrtdts*u_n.*dW1); % update u
    u_n=u_new;
end

localWN=zeros(J/S1+1,1);
for i=0:J/S1
    fftmatend=fft(u_new(J+2:end).*exp(-1/wl^2*(x+L-i*S1*h).^2)); %% 1D Fourier transform of every column of component j
    spectrum=abs((fftmatend(2:round(J/(2*S1)))));
    [~,localWN(i+1)]=max(spectrum(Start:End)); %% C contains maximum values per column, I contains rowindex where maximum is attained
    localWN(i+1)=localWN(i+1)+Start-1;
end
PulseCount=zeros(1,MAX);
PulseCount(Start:End)=histcounts(localWN,Start-0.5:1:End+0.5); %count how often each local wave number is found
end
