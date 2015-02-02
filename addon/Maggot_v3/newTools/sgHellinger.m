function hd = sgHellinger( f1, f2)

warning off MATLAB:divideByZero

m1=f1.mu;
s1=sqrt(f1.covariances);
m2=f2.mu;
s2=sqrt(f2.covariances);

hd=sqrt(1-sqrt(2*s1*s2/(s1^2+s2^2)) * exp(-(m1-m2)^2/4/(s1^2+s2^2)));
