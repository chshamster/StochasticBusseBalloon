function [ut]=spdeKlaus(IC,t,dt,J,a,sigma,q,dEE)
ut=zeros(2*(J+1),length(t));   %initialize vectors
ut(:,1)=IC; u_n=IC;            %set initial condition
sqrtdts=sqrt(dt)*sigma;

%compute settings for Lord's algorithms
tilde_q=[q; q(end-1:-1:2)];
Ntilde=length(tilde_q);
d=ifft(tilde_q,'symmetric')*Ntilde;
flag=false;
%time loop
for k=1:length(t)-1 
    if flag==false %generate two samples
      [dW11,dW21]=circ_cov_sampleB(d,Ntilde);
      dW1=cat(1,zeros(J+1,1),dW11(1:J+1)); dW2=cat(1,zeros(J+1,1),dW21(1:J+1));
      flag=true;
    else % use second sample from last call
      dW1=dW2; flag=false; 
    end
     fu=f(u_n,a,J); % evaluate f
     u_new=dEE\(u_n+dt*fu+sqrtdts*u_n.*dW1); % update u
     ut(:,k+1)=u_new; u_n=u_new; 
end
end
