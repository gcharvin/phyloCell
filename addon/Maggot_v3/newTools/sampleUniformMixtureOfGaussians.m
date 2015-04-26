function X = sampleUniformMixtureOfGaussians( centers, covariances, N )

lLim=min(centers-3*sqrt(covariances)');
uLim=max(centers+3*sqrt(covariances)');

X=lLim+rand(1,N)*(uLim-lLim);


