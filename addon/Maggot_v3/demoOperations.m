function demoOperations()

Nl = 30 ; 
dat = rand(3,Nl) ;
kde = executeOperatorIKDE( [], 'input_data', dat,'add_input' ) ;

x = rand(3,5) ;
selectSubDimensions = [] ;
res = executeOperatorIKDE( kde, 'input_data', x, 'evalPdfOnData', 'selectSubDimensions', selectSubDimensions ) ;
res.subRegularized
res.evalpdf

res = executeOperatorIKDE( kde, 'input_data', x, 'evalTypOnData', 'selectSubDimensions', selectSubDimensions ) ;
res.subRegularized
res.evaltyp
kde = res.kde ;

res = executeOperatorIKDE( kde, 'input_data', x, 'evalTypOnData', 'selectSubDimensions', selectSubDimensions ) ;

dat = rand(3,Nl) ;
selectSubDimensions = [] ;
kde2 = executeOperatorIKDE( [], 'input_data', dat,'add_input' ) ;
res = executeOperatorIKDE( kde, 'additional_kde', kde2, 'evalHellingerBetween' , 'selectSubDimensions', selectSubDimensions ) ;
res.subRegularized
res.distance_hell