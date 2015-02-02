function [mu1,sg1,n1]=updateMV(x,mu,sg,n,w)
%[mu1,sg1,n1]=updateMV(x,mu,sg,n)
%Update means and variances
%x: new data vector
%mu: current means
%sg: current variances
%n: current number of samples
%mu1: updated means
%sg1: updated variances
%n1: updated number of samples

if nargin<5
   w=1;
end

N=length(x);
mu1=zeros(1,N);
sg1=zeros(1,N);

if w>0
   for i=1:N
      % mu1(i)=1/(n+1)*(n*mu(i)+x(i));
      mu1(i) = 1/(n+w)*(n*mu(i) + w*x(i));
      % s=(n-1)*sg(i)+n*mu(i)*mu(i)+x(i)*x(i);
      s = (n-1)*sg(i) + n*mu(i)*mu(i) + w*x(i)*x(i);
      % sg1(i)=1/n*(s-(n+1)*mu1(i)*mu1(i));
      sg1(i) = 1/(n+w-1) * (s - (n+w)*mu1(i)*mu1(i));
   end;
   n1=n+w;
elseif w<0 && n>2 %&&  abs(x(1)-mu(1))<2*sqrt(sg(1))
   for i=1:N
         mu1(i)=1/(n-1)*(n*mu(i)-x(i));
         s=(n-1)*sg(i)+n*mu(i)*mu(i)-x(i)*x(i);
         sg1(i)=1/(n-2)*(s-(n-1)*mu1(i)*mu1(i));
         if sg1(i)<1e-9 
            mu1(i)=mu(i);
            sg1(i)=sg(i);
         end
   end;
   n1=n-1;
else
   mu1=mu;
   sg1=sg;
   n1=n;
end

