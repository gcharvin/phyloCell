function [input_kde, remaining_points] = addDatapointsByEMIfYouCan( input_kde, data, obs_mixing_weights, otherClasses )

diag_noise = 1e-5 ; %1e-10 ; %1e-8 ;% 1e-2 ; %1e-5 ;
pmin = 0.95 ; 0.95 ; % 0.99 ; 0.5 dela kul
% pmin2 = 0.5 ; 1e-4 ; 0.9 ; 
pmin2 = 0.7 ; 1e-4 ;
if size(data,2) > 1
    error('Designed only for one data-point!') ;
end

P = zeros(1,size(data,2)) ;
for i = 1 : length(otherClasses.pdfs)
    P = P + evaluatePointsUnderPdf(otherClasses.pdfs{i}, data,[], diag_noise)*otherClasses.priors*otherClasses.inner_priors(i) ;
end
prior_this = 1 - otherClasses.priors ;

% input_kde.pdf.smod.q(1).Mu = 'd' ; % to destroy smod!

[ p ] = evaluatePointsUnderPdf(input_kde.pdf, data, [],diag_noise)*prior_this ;
P = p./(p+P+1e-18) ;
id = P > pmin ;
remaining_points = data(:, id~=1) ;
points_to_insert = data(:, id) ;

if ~isempty(points_to_insert)
    
    [ idx, isvalid ] = getClosestComponent( input_kde.pdf, points_to_insert, pmin2, diag_noise ) ; %diag_noise  1e-1
%     isvalid = 1 ;
    if isvalid==0
        remaining_points = [remaining_points, points_to_insert] ;
        return ;
    end
    input_kde.pdf.w = input_kde.pdf.w*input_kde.ikdeParams.mix_weights(1) ;
    w_dat = obs_mixing_weights*input_kde.ikdeParams.mix_weights(2) ; 
 
        Mu = [input_kde.pdf.smod.q(idx).Mu, data] ;
        Cov = horzcat(input_kde.pdf.smod.q(idx).Cov, zeros(size(data,1))) ;         
       % w = [ input_kde.pdf.w(idx)*input_kde.pdf.smod.q(idx).w, w_dat ] ; 
       
       w = [ input_kde.pdf.smod.q(idx).w, w_dat/input_kde.pdf.w(idx) ] ; 
        
        [new_mu, new_Cov, w_out] = momentMatchPdf(Mu, Cov, w) ;
        input_kde.pdf.smod.q(idx).Mu = new_mu ;
        input_kde.pdf.q(idx).Cov = new_Cov ;
        input_kde.pdf.q(idx).w = w_out ; 
        input_kde.pdf.q(idx).w = input_kde.pdf.q(idx).w / sum(input_kde.pdf.q(idx).w) ;        

    input_kde.pdf.Mu(:,idx) = new_mu ;
    input_kde.pdf.Cov{idx} = new_Cov + input_kde.pdf.smod.H ;
    input_kde.pdf.w(idx) = w_out + input_kde.pdf.w(idx) ;    
    input_kde.pdf.smod.ps.Cov{idx} = new_Cov ; 
         
end
input_kde.pdf.w = input_kde.pdf.w / sum(input_kde.pdf.w) ;

% ---------------------------------------------------- %
function [ id, isvalid ] = getClosestComponent( mixtr, data, pminr, diag_noise )

Dn = eye(size(mixtr.Cov{1}))*diag_noise ;
dX = bsxfun(@minus, mixtr.Mu, data) ; 
p = zeros(1, length(mixtr.w)) ;
for i = 1 : length(mixtr.w)     
%     detD = det(mixtr.Cov{i}) ;
%     p(i) = (1/(sqrt(detD)))*exp(-0.5*sum(dX(:,i).*(mixtr.Cov{i}\dX(:,i)),1)) ; 
%     p(i) = (1/(sqrt(detD)))* exp(-0.5*dX(:,i)'*(mixtr.Cov{i}\dX(:,i))) ;
 p(i) =  exp(-0.5*dX(:,i)'*((mixtr.Cov{i}+Dn)\dX(:,i))) ;
%   p(i) = exp(-0.5*dX(:,i)'*((mixtr.Cov{i}+Dn)\dX(:,i))) ;
%  p(i) = exp(-0.5*dX(:,i)'*((mixtr.Cov{i}+Dn)\dX(:,i))) ;
end
[~, id] = max(p) ;% ptmp = p ; 
 isvalid = 0  ;
if max(p)  > pminr
    isvalid = 1 ;
 end

% p(id) = 1e-18 ;
% isvalid = 0  ;
% if max(p)/p(id) > pminr
%     isvalid = 1 ;
% end


% for i = 1 : length(mixtr.w)     
% %     detD = det(mixtr.Cov{i}) ;
% %     p(i) = (1/(sqrt(detD)))*exp(-0.5*sum(dX(:,i).*(mixtr.Cov{i}\dX(:,i)),1)) ; 
% %     p(i) = (1/(sqrt(detD)))* exp(-0.5*dX(:,i)'*(mixtr.Cov{i}\dX(:,i))) ;
% %  p(i) = (1/(sqrt(detD)))*exp(-0.5*dX(:,i)'*((mixtr.Cov{i}+Dn)\dX(:,i))) ;
%     p(i) = dX(:,i)'*((mixtr.Cov{i}+Dn)\dX(:,i)) ;
% %  p(i) = exp(-0.5*dX(:,i)'*((mixtr.Cov{i}+Dn)\dX(:,i))) ;
% end
% [~, id] = min(p) ;% ptmp = p ; 
%  isvalid = 0  ;
% if max(p) < 0.5
%     isvalid = 1 ;
%  end

