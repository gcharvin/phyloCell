function [actual_density,X]=my_marron_wand_normal_mixtures(density_number,N)

%
%    [actual_density,X]=
%            marron_wand_normal_mixtures(density_number,Y,N)
%     --------------------------------------------------------------------------
%     This function returns the actual density and samples from the
%     15 normal mixture densities used by Marron and Wand in
%      their paper:
%
%     J. S. Marron and M. P. Wand 'Exact Mean Integrated Squared
%     Error' The Annals of Statistics, 1992, Vol. 20, No. 2, 712-736
%     --------------------------------------------------------------------------
%     Author     :     Vikas.C.Raykar
%     Date        :     25 June 2005
%     Contact    :     vikas@cs.umd.edu
%     --------------------------------------------------------------------------
%     INPUTS
%     --------------------------------------------------------------------------
%     density_number --> can vary form 1-15 (See below)
%     Y                       --> data points at which the actual density is
%                                    to be evaluated
%     N                       --> number of samples.
%     --------------------------------------------------------------------------
%     OUTPUTS
%     --------------------------------------------------------------------------
%     actual_density   --> the actual density evaluated at Y.
%     X                      --->N samples from the density.
%     --------------------------------------------------------------------------
%
% 1   --> Gaussian
% 2   --> Skewed unimodal
% 3   --> Strongly skewed
% 4   --> Kurtotic unimodal
% 5   --> Outlier
% 6   --> Bimodal
% 7   --> Separated bimodal
% 8   --> Skewed bimodal
% 9   --> Trimodal
% 10 -->  Claw
% 11 -->  Double Claw
% 12 -->  Aymmetric Claw
% 13 -->  Asymmetric Double Claw
% 14 -->  Smooth Comb
% 15 -->  Discrete Comb

msg{1} = '1   --> Gaussian' ;
msg{2} = '2   --> Skewed unimodal' ;
msg{3} = '3   --> Strongly skewed' ;
msg{4} = '4   --> Kurtotic unimodal' ;
msg{5} = '5   --> Outlier' ;
msg{6} = '6   --> Bimodal' ;
msg{7} = '7   --> Separated bimodal' ;
msg{8} = '8   --> Skewed bimodal' ;
msg{9} = '9   --> Trimodal' ;
msg{10} = '10 -->  Claw' ;
msg{11} = '11 -->  Double Claw' ;
msg{12} = '12 -->  Aymmetric Claw' ;
msg{13} = '13 -->  Asymmetric Double Claw' ;
msg{14} = '14 -->  Smooth Comb' ;
msg{15} = '15 -->  Discrete Comb' ; 
% if ( density_number > 0 )
%     disp(sprintf('Sampling from: %s',msg{density_number})) ;
% end

Y = [] ;

if density_number==1
    weights = [1] ;
    m = [0] ;
    sigma = [1] ;
    actual_density=univariate_normal(0,1,Y);
    X=sample_univariate_normal(N,1,0,1,1);
end

if density_number==2
    G=3;
    weights=[1/5 1/5 3/5];
    m=[0 1/2 13/12];
    sigma=[1 2/3 5/9];
    actual_density=zeros(size(Y));
    for i=1:G
        actual_density=actual_density+weights(i)*univariate_normal(m(i),sigma(i),Y);
    end
    X=sample_univariate_normal(N,G,m,sigma,weights);
end

if density_number==3
    G=8;
    for j=1:8
        weights(j)=1/8;
        sigma(j)=(2/3)^(j-1);
        m(j)=3*((2/3)^(j-1)-1);
    end
    actual_density=zeros(size(Y));
    for i=1:G
        actual_density=actual_density+weights(i)*univariate_normal(m(i),sigma(i),Y);
    end       
    X=sample_univariate_normal(N,G,m,sigma,weights);
end

if density_number==4
    G=2;
    weights=[2/3 1/3];
    m=[0  0];
    mu = m ;
    sigma=[1 1/10];
    actual_density=zeros(size(Y));
    for i=1:G
        actual_density=actual_density+weights(i)*univariate_normal(m(i),sigma(i),Y);
    end
    X=sample_univariate_normal(N,G,m,sigma,weights);
end


if density_number==5
    G=2;
    weights=[1/10 9/10];
    m=[0  0];
    sigma=[1 1/10];
    actual_density=zeros(size(Y));
    for i=1:G
        actual_density=actual_density+weights(i)*univariate_normal(m(i),sigma(i),Y);
    end
    X=sample_univariate_normal(N,G,m,sigma,weights);
end

if density_number==6
    G=2;
    weights=[1/2 1/2];
    m=[-1  1];
    sigma=[2/3 2/3];
    actual_density=zeros(size(Y));
    for i=1:G
        actual_density=actual_density+weights(i)*univariate_normal(m(i),sigma(i),Y);
    end
    X=sample_univariate_normal(N,G,m,sigma,weights);
end

if density_number==7
    G=2;
    weights=[1/2 1/2];
    m=[-3/2  3/2];
    sigma=[2/3 2/3];
    actual_density=zeros(size(Y));
    for i=1:G
        actual_density=actual_density+weights(i)*univariate_normal(m(i),sigma(i),Y);
    end
    X=sample_univariate_normal(N,G,m,sigma,weights);
end

if density_number==8
    G=2;
    weights=[3/4 1/4];
    m=[0 3/2];
    sigma=[1 1/3];
    actual_density=zeros(size(Y));
    for i=1:G
        actual_density=actual_density+weights(i)*univariate_normal(m(i),sigma(i),Y);
    end
    X=sample_univariate_normal(N,G,m,sigma,weights);
end


if density_number==9
    G=3;
    weights=[9/20 9/20 1/10];
    m=[-6/5 6/5 0];
    sigma=[3/5 3/5 1/4];
    actual_density=zeros(size(Y));
    for i=1:G
        actual_density=actual_density+weights(i)*univariate_normal(m(i),sigma(i),Y);
    end
    X=sample_univariate_normal(N,G,m,sigma,weights);
end


if density_number==10
    G=6;
    weights(1)=1/2;
    sigma(1)=1;
    m(1)=0;
    for j=2:6
        weights(j)=1/10;
        sigma(j)=1/10;
        m(j)=((j-2)/2)-1;
    end
    actual_density=zeros(size(Y));
    for i=1:G
        actual_density=actual_density+weights(i)*univariate_normal(m(i),sigma(i),Y);
    end
    X=sample_univariate_normal(N,G,m,sigma,weights);
end

if density_number==11
    G=9;
    weights(1)=49/100;    sigma(1)=2/3;    m(1)=-1;
    weights(2)=49/100;    sigma(2)=2/3;    m(2)=1;
    for j=3:9
        weights(j)=1/350;
        sigma(j)=1/100;
        m(j)=((j-6)/2);
    end
    actual_density=zeros(size(Y));
    for i=1:G
        actual_density=actual_density+weights(i)*univariate_normal(m(i),sigma(i),Y);
    end
    X=sample_univariate_normal(N,G,m,sigma,weights);
end

if density_number==12
    G=6;
    weights(1)=1/2;
    sigma(1)=1;
    m(1)=0;
    k=2;
    for j=-2:1:2
        weights(k)=2^(1-j)/31;
        sigma(k)=2^(-j)/10;
        m(k)=j+(1/2);
        k=k+1;
    end
    actual_density=zeros(size(Y));
    for i=1:G
        actual_density=actual_density+weights(i)*univariate_normal(m(i),sigma(i),Y);
    end
    X=sample_univariate_normal(N,G,m,sigma,weights);
end


if density_number==13
    G=8;
    weights(1)=46/100;    sigma(1)=2/3;    m(1)=-1;
    weights(2)=46/100;    sigma(2)=2/3;    m(2)=1;
    k=3;
    for j=1:3
        weights(k)=1/300;
        sigma(k)=1/100;
        m(k)=-j/2;
        k=k+1;
    end
    for j=1:3
        weights(k)=7/300;
        sigma(k)=7/100;
        m(k)=j/2;
        k=k+1;
    end

    actual_density=zeros(size(Y));
    for i=1:G
        actual_density=actual_density+weights(i)*univariate_normal(m(i),sigma(i),Y);
    end
    X=sample_univariate_normal(N,G,m,sigma,weights);
end

if density_number==14
    G=6;
    k=1;
    for j=0:5
        weights(k)=(2^(5-j))/63;
        sigma(k)=((32/63))/(2^j);
        m(k)=(65-(96*((1/2)^(j))))/21;
        k=k+1;
    end

    actual_density=zeros(size(Y));
    for i=1:G
        actual_density=actual_density+weights(i)*univariate_normal(m(i),sigma(i),Y);
    end
    X=sample_univariate_normal(N,G,m,sigma,weights);
end


if density_number==15
    G=6;
    k=1;
    for j=0:2
        weights(k)=2/7;
        sigma(k)=2/7;
        m(k)=(12*j-15)/7;
        k=k+1;
    end
    for j=8:10
        weights(k)=1/21;
        sigma(k)=1/21;
        m(k)=(2*j)/7;
        k=k+1;
    end

    actual_density=zeros(size(Y));
    for i=1:G
        actual_density=actual_density+weights(i)*univariate_normal(m(i),sigma(i),Y);
    end
    X=sample_univariate_normal(N,G,m,sigma,weights);
end

if ( length(weights) == 1 )
    m = [m] ;
end
actual_density.norm.weights = weights ;
actual_density.norm.mu = m ;
actual_density.norm.covariances = (sigma').^2 ;


% figure;
% plot(Y,actual_density,'k');
% hold on;
% plot(X,0,'k+');

return

%--------------------------------------------------------------------------

function [density]=univariate_normal(m,sigma,Y)

density=(1/sqrt(2*pi*sigma*sigma))*exp(-((Y-m).^2)/(2*sigma*sigma));

return

%--------------------------------------------------------------------------

function [X]=sample_univariate_normal(N,G,m,sigma,weights)

% N_smp = zeros(1,length(weights)) ;
% cs = cumsum(weights) ;
% for i = 1 : N
%     f = cs <= rand(1) ; id = find(f) ;
%     if isempty(id) id = 1 ; else  id = id(length(id)) ; end
%     N_smp(id) = N_smp(id)+1 ;
% end

N_smp = zeros(1,length(weights)) ;
cs = cumsum(weights)*length(weights) ;
for i = 1 : N
    r =  rand(1)*length(weights)   ;
    f = r <= cs ; id = find(f) ;
    if isempty(id) id = 1 ; else  id = id(1) ; end
    N_smp(id) = N_smp(id)+1 ;
end



X = [] ;
for i = 1 : length(N_smp)
    if N_smp(i) > 0
        x = m(:,i)*ones(1,N_smp(i))+sigma(i)*randn(1,N_smp(i)); 
        X = [X,x] ;
    end
end
% 
% for i=1:G
%     n(i)=round(weights(i)*N);
% end
% n(G)=n(G)+(N-sum(n));
% 
% X=[];
% for i=1:G
%     X=cat(2,X,m(:,i)*ones(1,n(i))+sigma(i)*randn(1,n(i)));
% end

%--------------------------------------------------------------------------



