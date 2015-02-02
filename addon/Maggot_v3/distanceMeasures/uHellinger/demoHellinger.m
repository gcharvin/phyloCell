function demoHellinger()

% 
% HellError = 0.05 ;
% H_MC = []
% figure(1); clf ;
% N0 = 1 : 20 : 21 ;
% 
% for i = 1 : length(N0)
%     N = N0(i) ;
%     f0.mu = [  (rand(1,N)-0.5)*5]*0;
%     f0.weights=[  rand(1,N)/5 ] ;
%     f0.covariances = [  rand(1,N)*1.2 ]' ;
%     f0.weights = f0.weights/sum(f0.weights) ;
%     
%     Hmc = findOpt( f0, HellError ) 
%      
%     H_MC = [H_MC, Hmc] ;
%     figure(1); title(sprintf('Percent done: %f', 100*i/length(N0) )); drawnow ;
% end
% 
% plot(N0,H_MC) ;
% 
% function a = findOpt( f0, HellError ) 
% 
% f0.HellError = HellError ;
% a_initial = 1.5 ; % was 2.1 this is theoretical value for single gaussian at HellError=0.1
% options = optimset('Display','off','TolFun',0.001,'TolX',0.1,'LargeScale','on');
% % a = lsqnonlin(@a_function,a_initial,[1],[7],options,f0) ;
% a = fminbnd(@(a) a_function(a,f0),1,7,options) ;
% 
%  
%     fx = f0 ;
%     fx.covariances = fx.covariances*a ;
%     b = MCHellinger( fx, f0, 100000 ) ; 
%  
% if ( abs(b-f0.HellError) > 0.1 )
%     fsg = 67 ;
% end
% 
% function F = a_function(a,f0)
% 
% fx = f0 ;
% fx.covariances = fx.covariances*a ; 
% 
% F = 100*(f0.HellError - MCHellinger( fx, f0, 100000 ))^2 ; 
%  
%  

N = 6 ;
% f0.mu = [ 3.0991    0.6451   -3.2394    2.4057   -2.2278, (rand(1,N)-0.5)*5]*0   ;
% f0.weights=[0.1858    0.4807    0.0928    0.0020    0.2386, rand(1,N)/5 ] ; 
% f0.covariances = [0.947    1.844    0.0821    0.5006    0.143 , rand(1,N)*1.2+0.1 ]' ; % f0.covariances' = 0.6947    1.1844    0.5821    0.5006    0.7143

f0.mu = [(rand(1,N)-0.5)*5]*3   ;
f0.weights=[rand(1,N)/5 ] ; 
f0.covariances = [ rand(1,N)*1.2+0.1 ]' ;

% f0.mu = f0.mu(1) ;
% f0.covariances = f0.covariances(1) ;
% f0.weights = f0.weights(1);

f0.weights = f0.weights/sum(f0.weights) ;

 
%  f0.weights = ones(1, length(f0.weights)) ; 
%  f0.weights = f0.weights / sum(f0.weights) ;

% figure(1); clf ; 
% subplot(1,2,1) ; showDecomposedPdf(f0) ;

 H_MC = [] ; H_U = [] ;
 scales = [1.5 :0.1:3] ;
 
 scales = [3] ;
 for i = 1 : length(scales)

%      
%      scales(i) = 1.9 ;
 
R = [] ; B = [] ;


f0.mu = [(rand(1,N)-0.5)*5]*3   ;
f0.weights=[rand(1,N)/5 ] ; 
f0.covariances = [ rand(1,N)*1.2+0.1 ]' ;

% f0.mu = f0.mu(1) ;
% f0.covariances = f0.covariances(1) ;
% f0.weights = f0.weights(1);

f0.weights = f0.weights/sum(f0.weights) ;


     fx = f0 ; fx.covariances = f0.covariances*scales(i) ;
fx.mu = f0.mu + 3*rand(1,N) ;

for j = 1 :5 


    [ r_mc, r_ut] = MC_integrationVrtUtransform( fx, f0, 50000 ) ;
    Hu2 = uHellingerJointSupport2(fx, f0 ) ;
    Hu0 = uHellingerJointSupport(fx, f0 ) ;
    disp(sprintf('r_mc = %f, r_ut = %f, Hu_prejsnji = %f, Hu_zdaj = %f', r_mc, r_ut, Hu0, Hu2) ) ;
     
    disp(sprintf('Difference |r_mc-Hu_zdaj| = %f', abs(r_mc-Hu2) ) ) ;
%     disp(sprintf('Ratio Hu_prejsnji/Hu_zdaj = %f, Ratio |1-r_mc/Hu_zdaj| = %f',  Hu0/Hu2, abs(1 - r_mc/Hu2)  ) ) ;
    
    R = [R ,abs( r_mc - Hu2) ] ; 
    B = [B, Hu0/Hu2] ;
end
figure(1) ; clf ; plot(B);
[mean(R), max(R)]
[mean(B), max(B)]

%     Hu2/Hu0
return 
%      Hu = suHellinger(fx, f0 ) ; %
       Hu = uHellingerJointSupport2(fx, f0 ) ;
%        Hu =  sqrt( sqrt(2)*(Hu)  ) ;

%     Hu = getAvLikRatioAtSigmaPoints( fx, f0 )*2 ;
 
     N_inMonteCarlo = 10000 ;
     Hmc1 = calculateHellingerDistanceMonteCarlo( f0, fx, N_inMonteCarlo ) ;
     Hmc = MCHellinger( fx, f0, N_inMonteCarlo  ) ;
     
     [Hu , Hmc, Hmc1]
     
     H_U = [H_U, Hu] ; H_MC = [H_MC, Hmc] ;
 
     return ;
     figure(1); subplot(1,2,2) ; hold off ; plot(0,0,'.'); showDecomposedPdf(fx) ; drawnow ;
     
     figure(2); title(sprintf('Percent done: %f', 100*i/length(scales) )); drawnow ;
 end
 
 
 figure(3);  clf ; hold on ;
 plot(scales, H_MC, 'g') ; plot(scales,H_U, 'b') ;
% plot(scales,sqrt(2)*sqrt(1 - H_U), 'b') ;
 
 
%  
%  