function [ f ] = f(u,a,J)
f=zeros(size(u));
f(1:J+1)=a-u(1:J+1).*u(J+2:2*J+2).^2;   %set equation first component
f(J+2:2*J+2)=u(1:J+1).*u(J+2:2*J+2).^2; %set equation second component
end
