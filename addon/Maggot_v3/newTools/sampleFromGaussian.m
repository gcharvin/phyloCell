function r=sampleFromGaussian(n,mu,sg)
 
 %r=mu+sg*randn(1,n);
 
 sg=max(sg,1e-4);
 c=sqrt(3/2*sg);
 r=[mu-c*sg mu mu+c*sg];