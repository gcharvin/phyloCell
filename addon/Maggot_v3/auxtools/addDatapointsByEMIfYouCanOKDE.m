function [input_kde, remaining_points] = addDatapointsByEMIfYouCanOKDE( input_kde, data, rel_weights, obs_mixing_weights )

diag_noise = 1e-32 ; %1e-10 ; %1e-8 ;% 1e-2 ; %1e-5 ;
remaining_points = [] ;
% for i = 1 : size(data,2)
i = 1 ; % coded only for a single datapoint at a time! for more datapoints you'll have to rethink how
        % how to update the weights properly...
    [ idx, isvalid, dp ] = getClosestComponent( input_kde.pdf, data(:,i), diag_noise  ) ; %diag_noise  1e-1
    
    if isvalid        
        dst = bsxfun(@minus, input_kde.pdf.smod.q(idx).Mu, data(:,i)) ;
        dst = sum(dst.^2,1) ;
        [~, id_min] = min(dst) ;
    
            p2.Mu = [ input_kde.pdf.smod.q(idx).Mu(:,id_min), data(:,i)] ;
            p2.Cov = horzcat( input_kde.pdf.smod.q(idx).Cov(id_min), input_kde.pdf.smod.q(i).Cov{1}*0) ;
            p2.w =[rel_weights(1)*input_kde.pdf.w(idx)*input_kde.pdf.smod.q(idx).w(id_min), rel_weights(2)*obs_mixing_weights(i)] ;
                                        
            [new_mu, new_Cov, w_out] = momentMatchPdf(p2.Mu,p2.Cov, p2.w) ;
            input_kde.pdf.smod.q(idx).Mu(:,id_min) = new_mu ;
            input_kde.pdf.smod.q(idx).Cov{id_min} = new_Cov ;
            input_kde.pdf.smod.q(idx).w(id_min) = w_out ;
%               

            dmax = setdiff([1:length(input_kde.pdf.smod.q(idx).w)], id_min) ;
            if ~isempty(dmax)
               input_kde.pdf.smod.q(idx).w(dmax) = rel_weights(1)*input_kde.pdf.smod.q(idx).w(dmax)*input_kde.pdf.w(idx) ; 
            end

            input_kde.pdf.w = rel_weights(1)*input_kde.pdf.w ;
            input_kde.pdf.w(idx) = input_kde.pdf.w(idx) + rel_weights(2)*obs_mixing_weights(i) ;
            input_kde.pdf.smod.q(idx).w = input_kde.pdf.smod.q(idx).w / sum(input_kde.pdf.smod.q(idx).w) ;  
            
            % update upper layer
            [new_mu, new_Cov, w_out] = momentMatchPdf(input_kde.pdf.smod.q(idx).Mu,...
                                                      input_kde.pdf.smod.q(idx).Cov, ...
                                                      input_kde.pdf.smod.q(idx).w) ;
            input_kde.pdf.smod.ps.Cov{idx} = new_Cov ;
    else
        remaining_points = [remaining_points, data(:,i)] ;
    end    
% end
input_kde.pdf.w = input_kde.pdf.w / sum(input_kde.pdf.w) ;

% ---------------------------------------------------- %
function [ id, isvalid, dp ] = getClosestComponent( mixtr, data, diag_noise )

Dn = eye(size(mixtr.Cov{1}))*diag_noise ;
dX = bsxfun(@minus, mixtr.Mu, data) ; 
p = zeros(1, length(mixtr.w)) ;
% try
for i = 1 : length(mixtr.w)          
 p(i) = dX(:,i)'*((mixtr.Cov{i}+Dn)\dX(:,i)) ; 
end
% error(lastwarn) ;
% catch
%     df
% end
[dp, id] = min(p) ;% ptmp = p ; 
 isvalid = 0  ;
if  dp < 2.34  % max(p)  > pminr
    isvalid = 1 ;
end
 

