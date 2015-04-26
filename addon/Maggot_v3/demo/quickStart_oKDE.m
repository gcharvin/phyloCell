function quickStart_oKDE()
global pdfe
% 
% This is a skeleton you can use to get the oKDE [1] quickly up and running
% your own data.
%
% [1] Kristan Matej, Leonardis Ales and Skocaj Danijel, "Multivariate Online Kernel Density Estimation with Gaussian Kernels", 
% Pattern Recognition, 2011.
%
% Author: Matej Kristan, 2013.

% first install the oKDE toolbox
installEntireMaggot3() ; 

% we'll assume your data is in "dat", which is DxN matrix, D being the
% dimension of datapoints and N being the number of datapoints
D = 2 ; N = 100 ;
dat(1,:) = rand(1, N) ; % replace this with your own data!

dat(2,:) = rand(1, N)+4 ; % replace this with your own data!

% Note: if you have access to your datapoints in advance, or if you have a
% fair sample from your data points, you can use this to prescale your data
% prior to learning in the oKDE. This is especially convenient when the
% scale in one dimension (or subsspace) is significantly larger than in the
% other. Note that the oKDE will take care of this prescaling internally,
% but I still suggest that if you are able to provide some scaling factors
% in advance, you should do so.
prescaling = 1 ;
if prescaling
    [ Mu, T] = getDataScaleTransform( dat ) ;
    dat = applyDataScaleTransform( dat, Mu, T ) ;    
end

% initialize your KDE. Again, the oKDE has many valves to make it robust
% against poor initialization, but, if you can, it is recomended that you
% initialize it with sufficiently large sample set (N_init). A rule of thumb would
% be to initialize with more samples than twice the dimensionality of your data.

N_init = size(dat,1)*2 ; % how many samples will you use for initialization?
kde = executeOperatorIKDE( [], 'input_data', dat(:,1:N_init),'add_input' );

Dth = 0.02 ; % set the compression value (see the paper [1] for better idea about what value to use)
kde = executeOperatorIKDE( kde, 'compressionClusterThresh', Dth ) ;

% not you can add one sample at a time...
figure(1) ; clf ;
for i = N_init+1 : size(dat,2) 
    tic
    kde = executeOperatorIKDE( kde, 'input_data', dat(:,i), 'add_input'  ) ;
    t = toc ; 
    % print out some intermediate results
    msg = sprintf('Samples: %d ,Comps: %d , Last update time: %f ms\n', i , length(kde.pdf.w), t*1000 ) ; 
    title(msg) ; drawnow ;
end 
 
% your gaussian mixture model:
pdf_out = kde.pdf ;

% if you have prescaled your data, you will have to inverse the scaling of
% the estimated oKDE to project it back into the original space of your
% data. Note that from this point on, you can not update your pdf -- you
% will have to continue to update the KDE and inverse scaling again...
if prescaling    
  pdf_out = applyInvScaleTransformToPdf( pdf_out, Mu, T ) ;  
end

pdfe=pdf_out; 

% to evaluate a probability of a points x under the estimated kde:
%x = rand(D, 5) ; % some random data
%p = evaluatePointsUnderPdf(pdf_out, x ) 


% --------------------------------------------------------------------- %



