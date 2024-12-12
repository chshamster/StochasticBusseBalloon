clear all, close all

%This code computes solutions to the Klausmeijer model, Eq. (1.1) in
%Hamster, van Heijster and Siero, for sigma=0. An explanation of the numerical scheme can
%be found in Appendix A.2

%The folder Periodic Solutions contains a large set of stationary periodic
%solutions. 

run ../PeriodicSolutions/parameters.m   %Set parameters corresponding to the periodic solutions
a=1.45;     %rainfall parameters
k=30;       %wavenumber IC
load(['../PeriodicSolutions/n' num2str(k) 'a' num2str(a) 'profiles.mat']); %load periodic solution corresponding to a and k
profile(J+1,:)=profile(1,:); %set value at x=-L equal to x=L. 
IC=profile(:);               %The initial condition is the profile reshaped as vector


%Parameter settings for time stepping
N=2*10^4; T=50; dt=T/N; t=(0:dt:T)';

%Spatial discretization
Ab=computeAbP(J,h,d,1);             %Discretisation of the second derivative (in two blocks) with periodic BC
Lin=sparse(zeros(2*(J+1),2*(J+1))); %linear part of Klausmeijer model
Lin(1:J+1,1:J+1)=-speye(J+1,J+1);   
Lin(J+2:end,J+2:end)=-m*speye(J+1,J+1);
LinOp=Ab+Lin;                       %Full linear part of equation
EE=speye(2*(J+1))-dt*LinOp;         %Matrix to solve Eq. (A.4)
dEE=decomposition(EE);              %precompte decomposition to speed up performance

U=spdeKlausDet(IC,t,dt,J,a,dEE);    %solve equation

ut=U(1:J+1,:); vt=U(J+2:end,:);     %get both components
clear U;                            %Optional in case of memory issues

%speed up plotting
S1=16; S2=50;
utplot=ut(1:S1:end,1:S2:end);
vtplot=vt(1:S1:end,1:S2:end);

%plot solutions
figure(1)
subplot(1,2,1)
hold on
pcolor(x(1:S1:end),t(1:S2:end),utplot'); shading interp; colorbar;
title('$u$','Interpreter','Latex','Fontsize',30)
xlabel('$x$','Interpreter','Latex','Fontsize',40), xlim([-L L]);
ylabel('$t$','Interpreter','Latex','Fontsize',40)
hold off
subplot(1,2,2)
hold on
pcolor(x(1:S1:end),t(1:S2:end),vtplot'); shading interp; colorbar;
title('$v$','Interpreter','Latex','Fontsize',30)
xlabel('$x$','Interpreter','Latex','Fontsize',40), xlim([-L L]);
ylabel('$t$','Interpreter','Latex','Fontsize',40)
hold off

