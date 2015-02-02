%%
% Originally a part of: Maggot (developed within EU project CogX)
% Author: Matej Kristan, 2009 (matej.kristan@fri.uni-lj.si; http://vicos.fri.uni-lj.si/matejk/)
% Last revised: 2009
%%
function pdf1 = reestimateWeightsOfPdf( pdf1, pdf0 )

N_iterations = 50 ;

% pdf0 ... reference
e_min = 1e-5 ;
eps_dead = 1e-5 ;

% f_proposal = mergeDistributions( pdf1, pdf0, [0.5 0.5] ) ;

MaxV = 3 ;
[X, numSigPoints, w, k ] = getAllSigmaPointsOnMixture( pdf1, MaxV ) ;

sum(sum(imag(X)))

prior_weights = zeros(length(pdf1.w), size(X,2)) ;

% estimate a priori weights
for i = 1 : length(pdf1.w)
    prior_weights(i,:) = normpdf( X, pdf1.Mu(:,i), [], pdf1.Cov{i} ) ;
    prior_weights(i,:) =  prior_weights(i,:) / sum( prior_weights(i,:)) ;
    
%     idx = (i-1)*3+1 : (i-1)*3+3 ; prior_weights(i,:) = 0 ; prior_weights(i,idx) = 1 ; 
%     prior_weights(i,:) =  prior_weights(i,:) / sum( prior_weights(i,:)) ;
    
end
p_of_ref = evaluatePointsUnderPdf(pdf0, X) ;
p_of_ref0 = p_of_ref / sum(p_of_ref) ;

W1 = pdf1.w ;

for k = 1 : N_iterations
    for i = 1 : length(pdf1.w)
%         ss = 0 ;
%         for j = 1 : length(pdf0.w)
%           C = pdf1.Cov{i} + pdf0.Cov{j} ;
%           ss = ss + pdf0.w(j) * normpdf(pdf1.Mu(:,i), pdf0.Mu(:,j),[],C) ;
%         end
%         pdf1.w(i) = ss*sqrt(det(2*pdf1.Cov{i})) ;
%         pdf1.w = pdf1.w / sum(pdf1.w) ;
        
    [pdf1_off_k, pdf1_only_k] = splitK( pdf1, i ) ;
%        p_only_k = evaluatePointsUnderPdf(pdf1_only_k, X) ;
      p_only_k =  prior_weights(i,:) ;
%        p_only_k = p_only_k / sum(p_only_k) ;
       
        p_of_new = evaluatePointsUnderPdf(pdf1, X) ;
        idx = find(p_of_new < 1e-10) ;
        p_of_new(idx) = e_min ;
        W = (p_of_ref ./ p_of_new) .*p_only_k;
        w = sum(W)*pdf1.w(i) ;
        pdf1.w(i) = w ; 
        
        
%         pdf1.w
%         pdf1.w = pdf1.w / sum(pdf1.w) ;


%        wp = sum(p_only_k) ;
%        p_only_k = p_only_k / wp ; wp = 1 ;       
%        M = sum(X.*p_only_k)/sum(p_only_k) ; pdf1.Mu(:,i) ;
%    Mm = repmat(M,1,size(X,2)) ;
%       d = (Mm - X).*sqrt(p_only_k/sum(p_only_k)) ;
%        C = d*d'   ;
%        pdf1.Mu(:,i) = M ;
%        pdf1.Cov(i) = {C} ;
        
        pdf1.w = pdf1.w / sum(pdf1.w) ;
%         fignum = 1;
%     figure(fignum) ; clf ; 
%     drawDistributionGMM( 'pdf', pdf1, 'color', 'r' ) ;  
%     drawDistributionGMM( 'pdf', pdf0, 'color', 'g' ) ;
%     plot(X,X*0,'g*')
%     drawnow ; pause(0.1)
    end   
    pdf1 = removeDeadComponents( pdf1, eps_dead ) ;
    if length(pdf1.w) ~= length(W1)
        prior_weights = zeros(length(pdf1.w), size(X,2)) ;
        % estimate a priori weights
        for i = 1 : length(pdf1.w)
            prior_weights(i,:) = normpdf( X, pdf1.Mu(:,i), [], pdf1.Cov{i} ) ;
            prior_weights(i,:) =  prior_weights(i,:) / sum( prior_weights(i,:)) ;
        end
    else
        df_w = abs(W1 - pdf1.w) ;
        th = 1/length(W1) ;
        if ( max(df_w) < th*1e-3 )
            break ;
        end
    end
% [k max(df_w)]
    W1 = pdf1.w ;
    
%      pdf1.w = pdf1.w / sum(pdf1.w) ;
%  
%     fignum = 1;
%     figure(fignum) ; clf ; 
%     drawDistributionGMM( 'pdf', pdf1, 'color', 'r' ) ; 
%     drawDistributionGMM( 'pdf', pdf0, 'color', 'g' ) ;
%     plot(X,X*0,'g*'); 
%     msg = sprintf('Iteration: %d',k); title(msg)
%     drawnow ; %pause(0.5)
end


 

function [pdf1, pdf2] = splitK( pdf, k )

I = ones(1,length(pdf.w)) ;
I(k) = 0 ;
idx1 = find(I > 0) ; 
idx2 = k ;

pdf1.Mu = pdf.Mu(:,idx1) ;
pdf1.Cov = pdf.Cov(idx1) ;
pdf1.w = pdf.w(idx1) ;

pdf2.Mu = pdf.Mu(:,idx2) ;
pdf2.Cov = pdf.Cov(idx2) ;
pdf2.w = pdf.w(idx2) ;











return


% pdf0 ... reference
e_min = 1e-5 ;

f_proposal = mergeDistributions( pdf1, pdf0, [0.5 0.5] ) ;

MaxV = 3 ;
[X, numSigPoints, w, k ] = getAllSigmaPointsOnMixture( f_proposal, MaxV ) ;

sum(sum(imag(X)))

p_of_ref = evaluatePointsUnderPdf(pdf0, X) ;


prior_weights = zeros(length(pdf1.w), size(X,2)) ;
% estimate a priori weights
for i = 1 : length(pdf1.w)
    prior_weights(i,:) = normpdf( X, pdf1.Mu(:,i), [], pdf1.Cov{i} ) ; 
end

for it = 1 : 20
   for k = 1 : length(pdf1.w) 
       [pdf1_off_k, pdf1_only_k] = splitK( pdf1, k ) ;
       p_off_k = evaluatePointsUnderPdf(pdf1_off_k, X) ;
       p_only_k = evaluatePointsUnderPdf(pdf1_off_k, X) ;
                
       w = -sum((p_off_k - p_of_ref).*prior_weights(k,:)) / sum(p_only_k.^2) ;
       pdf1.w(k) = w ;
       pdf1.w = pdf1.w / sum(pdf1.w) ;
                   
       wp = sum(p_only_k) ;
       p_only_k = p_only_k / wp ; wp = 1 ;       
       M = sum(X.*p_only_k)/wp ;
       Mm = repmat(M,1,size(X,2)) ;
       d = (Mm - X).*sqrt(p/wp) ;
       C = d*d';
       pdf1.Mu(:,k) = M ;
       pdf1.Cov(k) = {C} ;
 
   end
   id = find(pdf1.w<0) ;
   pdf1.w(id) = 0 ;
   pdf1.w = pdf1.w / sum(pdf1.w) ;
   
   
   
   
   
   fignum = 1;
    figure(fignum) ; clf ; 
    drawDistributionGMM( 'pdf', pdf1, 'color', 'r' ) ; 
    drawDistributionGMM( 'pdf', pdf0, 'color', 'g' ) ;
    plot(X,X*0,'g*'); 
    msg = sprintf('Iteration: %d',it); title(msg)
    drawnow ; pause(0.5)
    
end




return ;
