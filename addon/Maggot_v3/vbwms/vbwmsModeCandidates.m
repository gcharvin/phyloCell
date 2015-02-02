%%
% Originally a part of: Maggot (developed within EU project CogX)
% Author: Matej Kristan, 2009 (matej.kristan@fri.uni-lj.si; http://vicos.fri.uni-lj.si/matejk/)
% Last revised: 2009
%%
function [converged, id_converged ] = vbwmsModeCandidates( data, weights, covariances, ...
                                                          precisions, determinants,...
                                                          stopThresh, search_area )
global plot_intermediate_results ;
  
num_datapoints = cols(data) ; 
dim = sqrt(cols(covariances)) ;
len_dim = dim^2 ;

Ctol = 1*eye(dim)*1e-15 ;

enforce_fixedBW = 0 ;
if ( ~isempty(search_area) ) 
    scale = 1; length(determinants) ;
    enforce_fixedBW = search_area.enforce_fixedBW ;
    Precision0 = getMinimalPrecision( precisions, determinants, scale ) ;
    if ~isempty(search_area.mu)
        enforce_search_area = 1 ;
    end
else
    enforce_search_area = 0 ;
end

if ( dim == 1 & plot_intermediate_results > 0 )
 D = max(data) - min(data) ;  
 bounds = [min(data)-D/3, max(data)+D/3] ;
 interval_pdf = [ bounds(1) : (bounds(2)-bounds(1))/100 : bounds(2) ] ; 
 prec_pdf = evaluateDistributionAt( data, weights, covariances, interval_pdf ) ;
 
 sign_p = 1 ; (-1)^(sum(prec_pdf) < 0) ;   
 
 prec_pdf_norm = 1 ; %1 / max(prec_pdf)
 prec_pdf = prec_pdf*prec_pdf_norm*sign_p ;
end
 
% generate kernel for stopping rule
a_weights = abs(weights) ; 
a_weights = a_weights / sum(a_weights) ;
% W = repmat(a_weights,len_dim,1)' ;
L_kern = sum(bsxfun(@times, precisions, a_weights'),1 ) ; %sum(precisions.*W,1) ;
% L_kern = bsxfun(@times, precisions, a_weights') ;
L_kern = reshape(L_kern,dim,dim) ;
L_kern = diag(L_kern) ;

% select all data points as initial centers
centers = data ;
% generate center index table
enumeration = [1:cols(centers)] ;
id_converged = [] ;

% precompute constants
constantsA = weights./sqrt(determinants) ;
%constantsA = abs(weights)/sum(abs(weights))./sqrt(determinants) ;
constantsB = data ;
PrecisionsR = {};
Covs = {} ;
for i_point = 1 : num_datapoints
    Precision = reshape(precisions(i_point,:),dim,dim) ;
    PrecisionsR = horzcat(PrecisionsR, Precision) ;
    constantsB(:,i_point) = Precision*data(:,i_point) ;    
end

% initialize converged centers
converged = [] ;

counter = 0 ;
 

while ~isempty( centers ) % continue while there are still some centers left
   prev_centers = centers ;
   num_centers = cols(centers) ;

%    W = zeros(num_centers, num_datapoints) ;
%    for i = 1 : num_datapoints
%    
%         D_2 = sqdist(centers,data(:,i),PrecisionsR{i})' ;
%  
%         W(:,i) = constantsA(i_point).* exp(-0.5*D_2) ;        
%    end
%    W = bsxfun(@times, W, 1./sum(W,1)) ;
%    
%    % update centers
%    for j = 1 : size(centers,2)  
%         centers(:,j) = inv(reshape(sum(bsxfun(@times, precisions, W(:,j) ),1),dim,dim)) * sum(bsxfun(@times, constantsB, W(j,:)),2) ;
%    end
%    
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   W = [] ;
   sumWeightPrec = zeros(num_centers,len_dim) ;
   sumB = centers*0 ;


   for i_point = 1 : num_datapoints
      point = data(:,i_point) ;
%       if enforce_fixedBW == 1 
%         Precision = reshape(Precision0,dim,dim) ;   
%       else
%        Precision = reshape(precisions(i_point,:),dim,dim) ;
%       end
%       Precision = reshape(precisions(i_point,:),dim,dim) ;
%       D_2 = sqdist(centers,point,Precision) ; 
        D_2 = sqdist(centers,point,PrecisionsR{i_point}) ; 
        W = constantsA(i_point).*exp(-0.5*D_2) ;
           
%       if enforce_fixedBW == 1 
%         Precision = repmat(Precision0,num_centers ,1 ) ;
%        else
%         Precision = repmat(precisions(i_point,:),num_centers ,1 ) ;
%       end
      
      % w*P^-1 term for all centers
%       Precision = repmat(precisions(i_point,:),num_centers ,1 ) ;
%       W1 = repmat(W',1,len_dim) ;  
       
%         sumWeightPrec = sumWeightPrec + W1.*Precision ; 
      

      % w*P^-1*x term for all centers
%       W2 = repmat(W,dim,1) ;
%       conB = repmat(constantsB(:,i_point),1,num_centers) ;
%       sumB = sumB + conB.*W2 ;
 
        sumWeightPrec = sumWeightPrec +  bsxfun(@times, precisions(i_point,:), W') ; 
        sumB = sumB + bsxfun(@times, constantsB(:,i_point), W) ;
   end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   counter = counter + 1 ;
 
   I0 = zeros(1, cols(centers)) ;
   % get new centers
   for i_center = 1 : cols(centers)
 
      centers(:,i_center) = inv(Ctol + reshape(sumWeightPrec(i_center,:),dim,dim))*sumB(:,i_center) ;
 
       
       
      % determine if center is out of search region  
      if enforce_search_area == 1
          dst = sqdist(search_area.mu, centers(:,i_center), search_area.inv_S) ;
          if ( dst > search_area.thresh )
            I0(i_center) = 1 ;
          end
      end
   end
   
   % determine converged centers to old centers
%    deltas = (prev_centers - centers).^2 ;
%    L_kernel = repmat(L_kern,1,cols(deltas)) ;
%    I = sqrt(sum(deltas.*L_kernel,1)) <= stopThresh ; 
    I = sqrt(sum(bsxfun(@times, (prev_centers - centers).^2, L_kern),1)) <= stopThresh ; 
 
   I = (I + I0)>0 ; 
      
   idx_converged = find(I) ;
   idx_non_converged = find(1-I) ;
   
   id_converged = [id_converged, enumeration(idx_converged)] ;
   enumeration = enumeration(idx_non_converged) ;
   converged = [converged,centers(:,idx_converged)] ;
   centers = centers(:,idx_non_converged) ;
 
   
 if plot_intermediate_results > 0  
    figure(plot_intermediate_results); clf ; 
    if ( dim == 1 ) 
        plot(interval_pdf,prec_pdf,'b') ; hold on ;
        pdf2 = evaluateDistributionAt( data, weights, covariances, centers ) ;
        pdf3 = evaluateDistributionAt( data, weights, covariances, converged ) ;
        plot(data,data*0,'*b') ;
        plot(centers, pdf2*prec_pdf_norm,'or') ;
        plot(converged,pdf3*prec_pdf_norm,'og') ;
        % hold on; plot(centers,centers*0+0.2,'or') ;
        % plot(converged,converged*0+0.2,'og') ;
        pause(0.1);
    elseif ( dim == 2 )
       plot(data(1,:),data(2,:),'*b') ; 
       hold on; plot(centers(1,:),centers(2,:),'or') ;
       plot(converged(1,:),converged(2,:),'og') ;       
    elseif ( dim == 3 )
       plot3(data(1,:),data(2,:),data(3,:),'*b') ; 
       hold on; plot3(centers(1,:),centers(1,:),centers(3,:),'or') ;
       plot3(converged(1,:),converged(2,:),converged(3,:),'og') ;     
    end
    msg = sprintf('Still active: %f',cols(centers))  ;
    title(msg) ; drawnow ;
 end
end

% ----------------------------------------------------------------------- %
function Precision0 = getMinimalPrecision( precisions, determinants, scale ) 

[val, id] = min(abs(determinants)) ;
Precision0 = precisions(id,:) * scale ;



