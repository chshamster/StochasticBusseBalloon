clear all, close all

%This code computes solutions to the Klausmeijer model, Eq. (1.1) in
%Hamster, van Heijster and Siero, and computes average exit times as function
%of wavenumber k and rainfall a. 
%Explanation of the numerical scheme can
%be found in Appendix A.2

%the folder Periodic Solutions contains a large set of stationary periodic
%solutions. 

%bunch of constants
run ../PeriodicSolutions/parametersc.m %set parameters of simulation
sigma=0.2;                             %set noise intensity

dt=0.02; Tmax=1e4;                   

ell=0.1; q=Cov(x,ell);                 %set correlation of the noise

%compute spatial discretization
Ab=computeAbP(J,h,d,1);                %Discretisation of second derivative (in two blocks)
Lin=sparse(zeros(2*(J+1),2*(J+1)));
Lin(1:J+1,1:J+1)=-speye(J+1,J+1);      %linear part of Klausmeijer model
Lin(J+2:end,J+2:end)=-m*speye(J+1,J+1);

LinOp=Ab+Lin;                          %Full linear part of equation
EE=speye(2*(J+1))-dt*LinOp;            %Matrix to solve Eq. (A.4)
clear Ab Lin LinOp                     %Optional for when menory issues arise

S2=200;                                %Every S2 timesteps, we compute the pulse number

K=25;                                  %Number of iterations for computing averages

%Discretization of the Busse balloon. For every a value, the matrix SE
%contains the lowest and highest wavenumber used. (SE --> Start-End).
Aa=0.40:0.05:2.85;
SE=[1,5;1,11;1,15;1,19;1,23;1,28;5,32;7,37;9,42;10,47;11,52;12,52;...
    13,52; 13,52; 14,52; 15,52; 15,52; 16,52; 17,51; 17,51; 18,51; 18,51; 19,51; 19,50; 19,50; 20,50; 20,50; 21,49; 21,49; 21,49; 21,49; 22,48;
    22,48; 22,48; 22,47; 22,47;22,47;22,46;22,46;22,45;22,45;23,44;23,44;24,43;25,42;26,42;27,41;28,40;30,39;31,37];
MAX=max(SE,[],'all'); 

dim=max(SE(:,2))-min(SE(:,1))+1; %number of different wavenumbers used per a value

%turn all indexes of a and SE values into a list to speed up the parrallel forloop.
List=[];
for j=1:length(Aa)
    for n=SE(j,1):SE(j,2)
        for i=1:K
            List=[List;[j,n,i]];
        end
    end
end

LL=length(List);
TexitList=zeros(LL,1);

%step through all combinations in List
parfor k=1:LL
    j=List(k,1);
    a=Aa(j);
    n=List(k,2);
    P=load (['../PeriodicSolutions/n' num2str(n) 'a' num2str(a) 'profiles.mat']);  %load initial condition
    profile=P.profile;
    profile(2*J+1,:)=profile(1,:);
    profile=interp1([-L:h/2:L],profile,x);  %interpolate initial condition to grid
    dEE=decomposition(EE); %precompute decomposition of EE for speed. (Cannot be done outside parfor) 
    TexitList(k)=spdeKlausIterTripleUntil(profile(:),dt,J,a,sigma,q,dEE,S2,n,Tmax);  %compute first exit time
end

%Turn list back into matrix
Texit=zeros(length(Aa),dim,K);
for k=1:LL
    Texit(List(k,1),List(k,2),List(k,3))=TexitList(k);
end

%Compute average Exit time
TexitM=sum(Texit,3)/K;

%Matrix containing pairs (k,a) that are used (useful for plotting)
InSim=nan(length(Aa),MAX);
for i=1:length(Aa)
    InSim(i,SE(i,1):SE(i,2))=1;
end

load ../plotm0/data.mat; %load boundary of Busse balloon
figure(1)
hold on
clim([0 Tmax])
myColorMap = parula;
colormap(myColorMap);
colorbar
fig2=image(Aa,1:MAX,(TexitM.*InSim)','CDataMapping','scaled');
set(fig2, 'AlphaData', ~isnan((TexitM.*InSim)'))
fimplicit(@(a,k) MostUnst(a,k,m,d,L),[0.9 2.88 0 100],'color','g','linewidth',3); %plot most unstable mode
plot(a1,500./WN1,'r','linewidth',3)
plot(a2,500./WN2,'r','linewidth',3)
plot(a3,500./WN3,'r','linewidth',3)
plot(a4,500./WN4,'r','linewidth',3)
plot(a5,500./WN5,'r','linewidth',3)
plot(a6,500./WN6,'r','linewidth',3)
xlabel('$a$','Interpreter','Latex','Fontsize',60)
ylabel('Wave Number I.C.','Interpreter','Latex','fontsize',60), ylim([0 MAX+1])
hold off


