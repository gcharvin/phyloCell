function demoMarginalization()
% examples of marginalizations

% learn a 3D GMM
Nl = 30 ; 
dat = rand(3,Nl) ;
kde = executeOperatorIKDE( [], 'input_data', dat,'add_input' ) ;

% evaluate data by marginalizing the second dimension
x = rand(3,5) ;
selectSubDimensions = [1 3] ;
res = executeOperatorIKDE( kde, 'input_data', x, 'evalPdfOnData', 'selectSubDimensions', selectSubDimensions ) ;
res.subRegularized
res.evalpdf

% marginalize out the second dimension and store the kde
selectSubDimensions = [1 3 ] ;
res = executeOperatorIKDE( kde, 'input_data', x, 'evalTypOnData', 'selectSubDimensions', selectSubDimensions ) ;
res.subRegularized
res.evaltyp
kde = res.kde ;

 