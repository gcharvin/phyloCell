%%
% Originally a part of: Maggot (developed within EU project CogX)
% Author: Matej Kristan, 2009 (matej.kristan@fri.uni-lj.si; http://vicos.fri.uni-lj.si/matejk/)
% Last revised: 2009
%%
function [model_new, ikdeParams] = manageAutomaticCompression( input_kde, otherClasses, obs )  
        


% now check if we should compress
% ikdeParams = input_kde.ikdeParams ;
%     if length(input_kde.pdf.w) > input_kde.ikdeParams.maxNumCompsBeforeCompression*1.3 
%         input_kde = executeOperatorIKDE( input_kde, 'compress_pdf' ) ;
%         ikdeParams.maxNumCompsBeforeCompression =  length(input_kde.pdf.w)*1.3 ; 
%     end
% model_new = input_kde.pdf ;
 
if input_kde.otherParams.MDL_memorylimitUseComps == 0
    
    if input_kde.otherParams.MDL_guides.MDL_guided_BW_Cprss == 1
        type = 'mdl' ;
    else
        type = 'histeresys' ;
    end
    
    switch type
        case 'mdl'
            scale_first = 1 ;
            scale_mid = 1 ;
            scale_last = 1 ;
        case 'histeresys'
            scale_first = 1.5 ;
            scale_mid = 0.5 ;
            scale_last = 0.6 ;
    end
    ikdeParams = input_kde.ikdeParams  ;
    
    N = abs(ikdeParams.maxNumCompsBeforeCompression) ;
    if length(input_kde.pdf.w) >  N - input_kde.ikdeParams.numComponentsAbsorbed  %|| rand(1) > 0.9 
        input_kde.ikdeParams.numComponentsAbsorbed = input_kde.ikdeParams.numComponentsAbsorbed/3 ;        
        
%         msg = sprintf('Compressing... num comps: %d', length(input_kde.pdf.w)) ; disp( msg ) ;
        input_kde = executeOperatorIKDE( input_kde, 'compress_pdf', 'otherClasses', otherClasses, 'input_data', obs ) ;
%         msg = sprintf('After compression.. num comps: %d', length(input_kde.pdf.w)) ; disp( msg ) ;
        len_new = length(input_kde.pdf.w) ;
%            disp('Compressing...')  
        if  ( len_new > N )
            N = N*scale_first ; %len_new*2 ;   %
        elseif len_new <= N*scale_mid
            N = N*scale_last ;
        end
        
        ikdeParams.maxNumCompsBeforeCompression = N ;
    end
    
else
   if length(input_kde.pdf.w) > input_kde.otherParams.MDL_memorylimitUseComps %input_kde.otherParams.MDL_memorylimit
        input_kde = executeOperatorIKDE( input_kde, 'compress_pdf', 'otherClasses', otherClasses ) ; 
   end
   ikdeParams = input_kde.ikdeParams  ;
end

model_new = input_kde.pdf ;
if isfield(input_kde, 'idxToref_out') && isfield(input_kde, 'addData')
   model_new.idxToref_out = input_kde.idxToref_out ;
   model_new.addData = input_kde.addData ;   
end



