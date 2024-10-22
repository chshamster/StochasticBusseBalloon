function [ut]=spdeKlausDet(IC,t,dt,J,a,dEE)
ut=zeros(2*(J+1),length(t)); % initialize vectors
ut(:,1)=IC; u_n=IC;          % set initial condition
for k=1:length(t)-1          % time loop
    fu=f(u_n,a,J);           % evaluate nonlinear part
    u_new=dEE\(u_n+dt*fu);   % solve for new time step
    ut(:,k+1)=u_new; u_n=u_new;
end
end