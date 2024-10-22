clear all, close all

%This code computes solutions to the Klausmeijer model, Eq. (1.1) in
%Hamster, van Heijster and Siero, and computes average exit times as function
%of wavenumber k and rainfall a. 
%Explanation of the numerical scheme can
%be found in Appendix A.2

%The folder Periodic Solutions contains a large set of stationary periodic
%solutions. 

%bunch of constants
run ../PeriodicSolutions/parametersc.m %set parameters of simulation
sigma=0.25;                             %set noise intensity

N=1e5; T=2500; dt=T/N; t=(0:dt:T)';                

ell=0.1; q=Cov(x,ell);                 %set correlation of the noise

%compute spatial discretization
Ab=computeAbP(J,h,d,1);                %Discretisation of second derivative (in two blocks)
Lin=sparse(zeros(2*(J+1),2*(J+1)));
Lin(1:J+1,1:J+1)=-speye(J+1,J+1);      %linear part of Klausmeijer model
Lin(J+2:end,J+2:end)=-m*speye(J+1,J+1);

LinOp=Ab+Lin;                          %Full linear part of equation
EE=speye(2*(J+1))-dt*LinOp;            %Matrix to solve Eq. (A.4)
clear Ab Lin LinOp                     %Optional for when menory issues arise

I=100;                                 %Number of iterations for computing averages

wl=50;                                 %window of the local wavenumber
S1=8;                                  %local wavenumber is computed for every S1 spatial steps

%Discretization of the Busse balloon. For every a value, the matrix SE
%contains the lowest and highest wavenumber used. (SE --> Start-End).
Aa=1.5:0.05:2.5;
SE=[ 18,34; 18,35; 19,36; 20,36; 20,36; 21,37; 21,37; 21,37; 22,37; 22,38;
   22,38; 22,38; 23,38; 23,38;23,38;23,38;24,38;24,38;24,38;24,39;25,39];

MAX=max(SE,[],'all');


%turn all indexes of a and SE values into a list to speed up the parrallel forloop.
List=[];
for i=1:length(Aa)
    for j=1:I
        List=[List; i,j];
    end
end
LL=length(List);

PulseCountList=zeros(LL,MAX);
parfor k=1:LL
    k
    a=Aa(List(k,1));
    %Start and End form the range in which the predominant local wavenumber
    %is found
    Start=SE(List(k,1),1);
    End=SE(List(k,1),2);
    dEE=decomposition(EE);
    FP=zeros(1,2); %compute stationary state
    FP(1)=(a/2-sqrt(a^2-4*m^2)/2); 
    FP(2)=(a-FP(1))/m;
    IC=[FP(1)*ones(size(x));FP(2)*ones(size(x))]; %set stationary state as IC
    PulseCountList(k,:)=spdeKlaus(IC,t,dt,L,x,J,h,a,sigma,q,dEE,S1,MAX,wl,Start,End);
end


%turn list back into matrix
PulseCount=zeros(length(Aa),I,MAX);
for i=1:LL
    PulseCount(List(i,1),List(i,2),:)=PulseCountList(i,:);
end

%compute average distribution
PulseCountAv=1/(J/S1+1)*squeeze(sum(PulseCount,2))/I;

%Matrix containing pairs (k,a) that are used (useful for plotting)
InSim=nan(length(Aa),MAX);
for i=1:length(Aa)
    InSim(i,SE(i,1):SE(i,2))=1;
end

load ../plotm0/data.mat;
figure(1)
hold on
myColorMap = parula;
colormap(myColorMap);
colorbar
fig=image(Aa,1:MAX,(PulseCountAv.*InSim)','CDataMapping','scaled');
set(fig, 'AlphaData', ~isnan((PulseCountAv.*InSim)'))
fimplicit(@(a,k) MostUnst(a,k,m,d,L),[1.475 2.525 0 100],'color','g','linewidth',3);
plot(a1,500./WN1,'r','linewidth',3)
plot(a2,500./WN2,'r','linewidth',3)
plot(a3,500./WN3,'r','linewidth',3)
plot(a4,500./WN4,'r','linewidth',3)
plot(a5,500./WN5,'r','linewidth',3)
plot(a6,500./WN6,'r','linewidth',3)
xlabel('$a$','Interpreter','Latex','Fontsize',50), xlim([1.4 2.9])
ylabel('Stationary Distribution','Interpreter','Latex','fontsize',50), ylim([0 MAX+1])
hold off
