%%
% Originally a part of: Maggot (developed within EU project CogX)
% Author: Matej Kristan, 2009 (matej.kristan@fri.uni-lj.si; http://vicos.fri.uni-lj.si/matejk/)
% Last revised: 2009
%%
function obs = sampleGaussianMixture( pdf, N, i )
% Draws N samples from a Gaussian mixture model pdf

if nargin > 2
    pdf.w = pdf.w(i) ; 
    pdf.w = pdf.w / sum(pdf.w) ;
    pdf.Mu = pdf.Mu(:,i) ;
    pdf.Cov = pdf.Cov(i) ;
end

if N < 10*length(pdf.w) 
%     N_smp = zeros(1,length(pdf.w)) ;
%     cs = cumsum(pdf.w) ;
%     for i = 1 : N
%         f = cs <= rand(1) ; id = find(f) ;
%         if isempty(id) id = 1 ; else  id = id(length(id)) ; end 
%         N_smp(id) = N_smp(id)+1 ;
%     end
    N_smp = zeros(1,length(pdf.w)) ;
    cs = cumsum(pdf.w)*length(pdf.w) ;
    for i = 1 : N
        r =  rand(1)*length(pdf.w)   ;
        f = r <= cs ; id = find(f) ;
        if isempty(id) id = 1 ; else  id = id(1) ; end 
        N_smp(id) = N_smp(id)+1 ;
    end
    N0 = N ;
else
    N0 = N ;
    N_smp = round((pdf.w)*N) ;
end
% N = 100 ;
% N_smp = round((pdf.w)*N) ;
% [a,a_i] = max(N_smp) ;
% N_smp(a_i) = N_smp(a_i) - (sum(N_smp)-N) ;

x = [] ;
obs = [] ;
% I = [] ;
for i = 1 : length(pdf.w)
    if (N_smp(i) > 0)        
%         I = [I,i]
        x = gaussSample(pdf.Mu(:,i), pdf.Cov{i}, N_smp(i))' ;
    end
    obs = [obs, x]  ;
end
obs = obs(:,randperm(min(N0,size(obs,2)))) ;
% ----------------------------------------------------------------------- %
function x = gaussSample(mu, covar, nsamp)
%GSAMP	Sample from a Gaussian distribution.
%
%	Description
%
%	X = GSAMP(MU, COVAR, NSAMP) generates a sample of size NSAMP from a
%	D-dimensional Gaussian distribution. The Gaussian density has mean
%	vector MU and covariance matrix COVAR, and the matrix X has NSAMP
%	rows in which each row represents a D-dimensional sample vector.
%
%	See also
%	GAUSS, DEMGAUSS
%

%	Copyright (c) Ian T Nabney (1996-2001)

d = size(covar, 1);

mu = reshape(mu, 1, d);   % Ensure that mu is a row vector

% [evec, eval] = eig(covar);
[evec,eval,v]=svd(covar) ;

coeffs = randn(nsamp, d)*sqrt(eval);

x = ones(nsamp, 1)*mu + coeffs*evec';