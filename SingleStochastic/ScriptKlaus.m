clear all, close all
%This code computes a single realization to the Klausmeijer model, Eq. (1.1) in
%Hamster, van Heijster and Siero. 
%Explanation of the numerical scheme can
%be found in Appendix A.2

%The folder Periodic Solutions contains a large set of stationary periodic
%solutions. 

%bunch of constants
run ../PeriodicSolutions/parameters.m  
a=1.45;     %rainfall parameters
k=30;       %wavenumber IC
load(['../PeriodicSolutions/n' num2str(k) 'a' num2str(a) 'profiles.mat']); %load periodic solution corresponding to a and k
profile(J+1,:)=profile(1,:); %set value at x=-L equal to x=L. 
IC=profile(:);               %The initial condition is the profile reshaped as vector

sigma=0.05; %noise intentisy

%settings for discretizations
N=6*10^4; T=500; dt=T/N; t=(0:dt:T)';

ell=0.1;       %correlation length
q=Cov(x,ell);  %correlation function
Ab=computeAbP(J,h,d,1);             %Discretisation of the second derivative (in two blocks) with periodic BC
Lin=sparse(zeros(2*(J+1),2*(J+1))); %linear part of Klausmeijer model
Lin(1:J+1,1:J+1)=-speye(J+1,J+1);   
Lin(J+2:end,J+2:end)=-m*speye(J+1,J+1);
LinOp=Ab+Lin;                       %Full linear part of equation
EE=speye(2*(J+1))-dt*LinOp;         %Matrix to solve Eq. (A.4)
dEE=decomposition(EE);              %precompte decomposition to speed up performance

%Solve equation
U=spdeKlaus(IC,t,dt,J,a,sigma,q,dEE);
ut=U(1:J+1,:); vt=U(J+2:end,:); %get both components
clear U;                        %usefull when storing the solution is an issue    

%speed up plotting
S1=16; S2=10;
utplot=ut(1:S1:end,1:S2:end);
vtplot=vt(1:S1:end,1:S2:end);
clear ut vt                     %optional when storing the solution is an issue


figure(1)
hold on
pcolor(x(1:S1:end),t(1:S2:end),utplot'); shading interp; colorbar;
xlabel('x','Interpreter','Latex','Fontsize',40), xlim([-L L]);
ylabel('t','Interpreter','Latex','Fontsize',40)
hold off


figure(2)
hold on
pcolor(x(1:S1:end),t(1:S2:end),vtplot'); shading interp; colorbar;
xlabel('x','Interpreter','Latex','Fontsize',40), xlim([-L L]);
ylabel('t','Interpreter','Latex','Fontsize',40)
hold off


