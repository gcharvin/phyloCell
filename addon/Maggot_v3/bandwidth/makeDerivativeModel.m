%%
% Originally a part of: Maggot (developed within EU project CogX)
% Author: Matej Kristan, 2009 (matej.kristan@fri.uni-lj.si; http://vicos.fri.uni-lj.si/matejk/)
% Last revised: 2009
%%
function derivativeModel = makeDerivativeModel( model, obs, varargin  )

SourceBw = NaN ;
obs_mixing_weights = NaN ;
mix_weights = NaN ;
N_eff = NaN ;
covNew = NaN ;
typeDerivative = 'normalRule' ;
% process arguments
args = varargin;
nargs = length(args);
for i = 1:2:nargs
    switch args{i}
        case 'SourceBw', SourceBw = args{i+1} ;
        case 'obs_mixing_weights', obs_mixing_weights = args{i+1} ;
        case 'mix_weights', mix_weights = args{i+1} ;
        case 'typeDerivative', typeDerivative = args{i+1} ;
        case 'covNew', covNew = args{i+1} ;
        case 'N_eff', N_eff = args{i+1} ;
        otherwise 
            msg = sprintf('Unknown parameter: %s !', args{i}) ;
            error(msg) ;
    end
end

% verify if there are any components left in the model and exit if so.
if ( length(model.w) == 0 )
    derivativeModel = [] ;
    return ;
end


if isequal(typeDerivative,'normalRule')
    % compute covariance using the Scot's rule of thumb
    d = size(model.Cov{1},1) ;
    H_est = covNew*(N_eff)^(-2/(d+4)) ;
    if isnan(covNew)
        error('Data covariance matrix not given!') ;
    end
elseif isequal(typeDerivative,'useSourceBw')
    H_est = SourceBw ;
else
   error('Unknown value for parameter typeDerivative!!!!') ;
end

derivativeModel = augmentMixtureModelByThis( model, obs, H_est, obs_mixing_weights,...
                                             mix_weights ) ;