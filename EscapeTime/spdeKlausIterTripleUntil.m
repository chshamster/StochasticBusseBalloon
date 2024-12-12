function [Texit]=spdeKlausIterTripleUntil(IC,dt,J,a,sigma,q,dEE,S2,n,Tmax)
u_old=IC; % set initial condition
sqrtdts=sqrt(dt)*sigma;
%compute settings for Lord's algorithms
tilde_q=[q; q(end-1:-1:2)];
Ntilde=length(tilde_q);
d=ifft(tilde_q,'symmetric')*Ntilde;
flag=false;
k=2;
pn=n; %set pulse number to the number of the IC
T=0;
while pn==n && T<Tmax    %time loop. Keep running until the pulse number changes. 
    if flag==false %generate two samples
        [dW11,dW21]=circ_cov_sampleB(d,Ntilde);
        dW1=[zeros(J+1,1);dW11(1:J+1)]; dW2=[zeros(J+1,1);dW21(1:J+1)];
        flag=true;
    else %use second sample from last call
        dW1=dW2; flag=false;
    end
    fu=f(u_old,a,J); % evaluate f
    u_new=dEE\(u_old+dt*fu+sqrtdts*u_old.*dW1); %update u
    if mod(k-1,S2)==0  %compute pulse number every S2 steps
        smoothvt=smoothdata(u_new(J+2:2*J+2),'gaussian',32);
        pn=max(sum(islocalmax(smoothvt,'MinProminence',0.3)),sum(islocalmin(smoothvt,'MinProminence',0.3)));
    end
    u_old=u_new;
    k=k+1;
    T=T+dt;
end
Texit=T;
end



