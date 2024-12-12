function [Ab] = computeAbP(J,h,dF,dR)
%discretization of the second derivative in two blocks with periodic boundary
%conditions
e = ones(J+1,1); A = spdiags([-e 2*e -e], -1:1, J+1, J+1);
A(1,end)=-1; A(end,1)=-1;
A=-1/h^2*A;
Ab=zeros(2*(J+1),2*(J+1));
Ab(1:J+1,1:J+1)=dF*A; Ab(J+2:2*J+2,J+2:2*J+2)=dR*A; 
Ab=sparse(Ab);
end

