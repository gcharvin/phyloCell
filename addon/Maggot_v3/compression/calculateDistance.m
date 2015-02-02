%%
% Originally a part of: Maggot (developed within EU project CogX)
% Author: Matej Kristan, 2009 (matej.kristan@fri.uni-lj.si; http://vicos.fri.uni-lj.si/matejk/)
% Last revised: 2009
%%
function [d , costThreshold] = ...
     calculateDistance(pdf, costFunction, costThreshold,...
                       numberOfSamples, MDL_params, useMargHellingerCompression)


switch costFunction
    case 'hellinger'
        d = distanceHellingerPdfs( pdf, useMargHellingerCompression ) ; 
    case 'numberOfComponents'
        d = length(pdf.w) < costThreshold ;
    case 'alpha_MDL'
         MDL = MDLbetweenDistributions( 'input_params', MDL_params, ...
                                        'pdf', pdf ) ;         
        if isempty(costThreshold)
           costThreshold = MDL ; 
        end
                 
        if MDL < costThreshold % accept distribution
            costThreshold = MDL ; 
            d = costThreshold - 1 ;        
        else % reject this solution
            d = costThreshold + 1 ;
        end
end

% ------------------------------------------------------------------ %
function d = distanceHellingerPdfs( pdf, useMargHellingerCompression )
global reference_pdf ;
 
d = uHellingerJointSupport2_ND( pdf, reference_pdf,...
                                'useMarginals', useMargHellingerCompression ) ;

% figure(2); clf; subplot(1,2,1); drawDistributionGMM( 'pdf', pdf2 ) ; subplot(1,2,2); drawDistributionGMM( 'pdf', reference_pdf ) ;
% msg = sprintf('Distance: %1.5g',d) ; title(msg) ; 
% pause
% 
%------------------------------------------------------------------------ %