
function H = uCostModel_full( negModel, posModel, posModel_r, negModelPrior, use_approximate_calc, pdf_for_sigmas, type_cost )
% Matej Kristan (2009)
% calculates a cost of model reduction in terms of classification
% accuracy.
 

% od testa za Letter se ta razlikuje po tem, da sem aktiviral brisanje
% predhranjenja in doloèil v funkciji ki klièe to, da naj raèuna na toèkah
% s subseta, ki ga želimo kompresirati namesto na celotni pozitivni pdf.

% focus only on centers
% focus_on_centers = 0 ; % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

% use_approximate_calc = 1 ;
use_sigma_points = 0 ;

% if use_sigma_points == 1 
%     focus_on_centers = 0 ;
% else
%     focus_on_centers = 1 ;
% end

% check if negative model even exists for the local compression:
emptynegative = 0 ;
if isempty(negModel)  
    emptynegative = 1 ;
end
if emptynegative == 1 
   
    if isempty(posModel_r)
        H.sigmaPoints = [] ;
        H.precalcs = [] ;
    else
         H = 0 ;
    end
    return ;
end 

% negModel.precalcStat = [] ;

modelPriors.pNeg = negModelPrior ; % 0.5 ;
modelPriors.pPos = 1 - modelPriors.pNeg ; 

modelPriors.pPos = 0.5 ;
modelPriors.pNeg = 0.5 ;
% kasneje še negModel.inner_priors(i)=0.5 !!!!!!!!!!!!!!!!!!!!!!!!!

minTol = 1e-50 ;
nimPerc = 0.001 ;

% for i = 1 : length(negModel.pdfs)
%         negModel.inner_priors(i) = 0.5 ;
% end

modelPriors.pPos = 1 / (length(negModel.pdfs) + 1) ;
modelPriors.pNeg = 1 - modelPriors.pPos ; 
 
if ~isfield(negModel,'precalcStat') || isempty(negModel.precalcStat)
   % generate sigma points 
    MaxV = 3 ; 
    f0 = posModel ;
    
    % calculate sigma points for entire set
    if use_sigma_points == 1
       [X, sigPointsPerComponent, w, k ] = getAllSigmaPointsOnMixture( f0, MaxV ) ;        
        X = real(X) ;
        W = repmat(f0.w,sigPointsPerComponent,1) ;
        W = reshape(W,1,length(f0.w)*sigPointsPerComponent) ;
        w2 = repmat(w,1,length(f0.w)) ;
        W = W.*w2 ;
    else
%         X = f0.Mu ;
%         W = f0.w ; 
    end
    % store vectors
    sigmaPoints.X = [f0.Mu] ;
    sigmaPoints.W = [f0.w] ;
    n_pdf.Mu = [] ; n_pdf.Cov = {} ; n_pdf.w = [] ;
    for i = 1 : length(negModel.pdfs)    
%         sigmaPoints.X = [sigmaPoints.X,negModel.pdfs{i}.Mu ] ;
%         sigmaPoints.W = [sigmaPoints.W, negModel.pdfs{i}.w] ;  
        n_pdf.Mu = horzcat(n_pdf.Mu, negModel.pdfs{i}.Mu) ;
        n_pdf.Cov = horzcat(n_pdf.Cov, negModel.pdfs{i}.Cov) ;
        n_pdf.w = horzcat(n_pdf.w, negModel.pdfs{i}.w) ;
        negModel.inner_priors(i) = 1/(length(negModel.pdfs)+1) ; 
    end
    n_pdf.w = n_pdf.w / sum(n_pdf.w) ;
    sigmaPoints.W = sigmaPoints.W / sum(sigmaPoints.W) ;
    ppos = evaluatePointsUnderPdf(posModel, sigmaPoints.X ) ;
    modelPriors.pPos = 1 / (length(negModel.pdfs) + 1) ;
    negModel.inner_priors = length(negModel.pdfs) / (length(negModel.pdfs) + 1) ;
    
    p1_z_c0 = ppos*modelPriors.pPos ;
    p1_z_c1 = evaluatePointsUnderPdf(n_pdf, sigmaPoints.X)*negModel.inner_priors ;
    p1_z = p1_z_c0 + p1_z_c1;
    p1_c0_g_z = p1_z_c0 ./ (p1_z + minTol) ;
    
    precalcs.p1_z_c1 = p1_z_c1 ; % negative model evaluated
    precalcs.ppos = ppos ; % positive referece model evaluated posterior
    precalcs.p1_c0_g_z = p1_c0_g_z ;
    
    H.sigmaPoints = sigmaPoints ;
    H.precalcs = precalcs ;
    if isempty(posModel_r)
        return ;
    else
        negModel.precalcStat = H ;
    end
end
 
if isempty(negModel.precalcStat.sigmaPoints)
    H = 0 ;
    return ;
end

X = negModel.precalcStat.sigmaPoints.X ;
W = negModel.precalcStat.sigmaPoints.W ;
p2_z_c0 = evaluatePointsUnderPdf(posModel_r, X)*modelPriors.pPos ;
p1_c0_g_z = negModel.precalcStat.precalcs.p1_c0_g_z ;
p1_c1_g_z = 1 - p1_c0_g_z ;
p2_z_c1 = negModel.precalcStat.precalcs.p1_z_c1 ;
p2_z = p2_z_c0 + p2_z_c1 ;
p2_c0_g_z = p2_z_c0 ./ (p2_z + minTol) ;
p2_c1_g_z = 1 - p2_c0_g_z ; %p2_z_c1 ./ (p2_z + minTol) ;
    
g0 = ((sqrt(p1_c0_g_z) - sqrt(p2_c0_g_z)).^2) ;
g1 = ((sqrt(p1_c1_g_z) - sqrt(p2_c1_g_z)).^2) ;
Hc = sum(W.*(g0 + g1)) ;
H = sqrt(Hc/2) * length(negModel.pdfs) ;

% H = sum(W.*sqrt((g0 + g1)/2)) ;
% H = sqrt(Hc/2) * length(negModel.pdfs) ;

% H = max(H_tmp) * length(negModel.precalcStat.precalcs.p1_z_c1) ;
if isempty(H)
    H = 0 ;
end

% ----------------------------------------------------------------------- %
function H = uCostModel_x( negModel, posModel, posModel_r, negModelPrior, use_approximate_calc, pdf_for_sigmas, type_cost )
% ta je tista, kjer je H=0.1/Nloc !!
% Matej Kristan (2009)
% calculates a cost of model reduction in terms of classification
% accuracy.
 

% od testa za Letter se ta razlikuje po tem, da sem aktiviral brisanje
% predhranjenja in doloèil v funkciji ki klièe to, da naj raèuna na toèkah
% s subseta, ki ga želimo kompresirati namesto na celotni pozitivni pdf.

% focus only on centers
% focus_on_centers = 0 ; % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

% use_approximate_calc = 1 ;
use_sigma_points = 0 ;
use_sampler = 0 ;
type_cost = 1 ;

% if use_sigma_points == 1 
%     focus_on_centers = 0 ;
% else
%     focus_on_centers = 1 ;
% end

% check if negative model even exists for the local compression:
emptynegative = 0 ;
if isempty(negModel)  
    emptynegative = 1 ;
end
if emptynegative == 1 
   
    if isempty(posModel_r)
        H.sigmaPoints = [] ;
        H.precalcs = [] ;
    else
         H = 0 ;
    end
    return ;
end 

% negModel.precalcStat = [] ;

modelPriors.pNeg = negModelPrior ; % 0.5 ;
modelPriors.pPos = 1 - modelPriors.pNeg ; 

modelPriors.pPos = 0.5 ;
modelPriors.pNeg = 0.5 ;
% kasneje še negModel.inner_priors(i)=0.5 !!!!!!!!!!!!!!!!!!!!!!!!!

minTol = 1e-50 ;
nimPerc = 0.001 ;


for i = 1 : length(negModel.pdfs)
        negModel.inner_priors(i) = 0.5 ;
end
modelPriors.pPos = 0.5 ;
modelPriors.pNeg = 0.5 ;

% nneg =length(negModel.pdfs) ;
% for i = 1 : length(negModel.pdfs)
%         negModel.inner_priors(i) = 1/(nneg+1) ;
% end
% modelPriors.pPos = nneg/(nneg+1) ;
% modelPriors.pNeg = 1/(nneg+1) ;
 
if ~isfield(negModel,'precalcStat') || isempty(negModel.precalcStat) ||1==1
   % generate sigma points 
    MaxV = 3 ; 
    f0 = posModel ;
    
    % calculate sigma points for entire set
    if use_sigma_points == 1
       [X, sigPointsPerComponent, w, k ] = getAllSigmaPointsOnMixture( f0, MaxV ) ;        
        X = real(X) ;
        W = repmat(f0.w,sigPointsPerComponent,1) ;
        W = reshape(W,1,length(f0.w)*sigPointsPerComponent) ;
        w2 = repmat(w,1,length(f0.w)) ;
        W = W.*w2 ;
    else
        
        if use_sampler ~= 1
%             if isempty(posModel_r)
                X = [f0.Mu  ];
                W = [f0.w  ] ; 
%             else
%                 X = [f0.Mu , posModel_r.Mu];
%                 W = [f0.w , posModel_r.w ]*0.5; 
%             end
        else
            fx = f0 ;
             for i = 1 : length(negModel.pdfs) 
                 fx.Mu = horzcat(fx.Mu, negModel.pdfs{i}.Mu) ;
                 fx.Cov = horzcat(fx.Cov, negModel.pdfs{i}.Cov) ;
                 fx.w = horzcat(fx.w, negModel.pdfs{i}.w) ;
             end
             fx.w  = fx.w  / sum(fx.w) ;
            
             [new_mu, new_Cov, w_out] = momentMatchPdf(fx.Mu, fx.Cov, fx.w) ;
             fx.Mu = new_mu ;
             fx.Cov = {new_Cov} ;
             fx.w = 1 ;
             [X, sigPointsPerComponent, w, k ] = getAllSigmaPointsOnMixture( fx, MaxV ) ;
             
             X = real(X) ;
             W = repmat(fx.w,sigPointsPerComponent,1) ;
             W = reshape(W,1,length(fx.w)*sigPointsPerComponent) ;
             w2 = repmat(w,1,length(fx.w)) ;
             W = W.*w2 ;
            %
%              x = gaussSample(new_mu, new_Cov, size(f0.Mu,2)*size(f0.Mu,1)*2) ;
%             
%             X = sampleGaussianMixture( f0, size(f0.Mu,2)*size(f0.Mu,1)*2  ) ;
%             W = ones(1,size(X,2),1)/size(X,2) ;
        end
    end
    % store vectors
    sigmaPoints.X = {} ;
    sigmaPoints.W = {} ;

    %ppos = evaluatePointsUnderPdf(posModel, X) ;
       
    % evaluate probabilities of model 1
    %p1_z_c0 = ppos*modelPriors.pPos ; 
 
    p1_z_c1 = {} ; p1_z = {} ; p1_c0_g_z = {} ; %p1_c1_g_z = {} ;
    ppos = {} ; p1_z_c0 = {} ;
    for i = 1 : length(negModel.pdfs)   
        if use_sampler ~= 1
%             X_tmp = [X,negModel.pdfs{i}.Mu] ;
%             W_tmp = [W,negModel.pdfs{i}.w ]*0.5 ;
            nneg =length(negModel.pdfs) ;
            X_tmp = [X,negModel.pdfs{i}.Mu] ;
            W_tmp = [W*modelPriors.pPos,negModel.pdfs{i}.w*modelPriors.pNeg ] ;
            W_tmp = W_tmp / sum(W_tmp) ;
%                 X_tmp = [X] ;
%                 W_tmp = [W] ;
        else
%             Xt = sampleGaussianMixture( f0, size(f0.Mu,2)*size(f0.Mu,1)*2  ) ;
%             Wt = ones(1,size(Xt,2),1)/size(Xt,2) ;
%             
%             X_tmp = [X,Xt] ;
%             W_tmp = [W,Wt]*0.5 ;
                X_tmp = [X] ;
                W_tmp = [W] ;
        end
        
        
        
        ppos_tmp = evaluatePointsUnderPdf(posModel, X_tmp) ;
        
        ppos = horzcat(ppos, ppos_tmp) ;
        p1_z_c0 = horzcat(p1_z_c0, ppos_tmp*modelPriors.pPos) ;
        
        negModel.inner_priors(i) = 0.5 ;        
        p_tmp = evaluatePointsUnderPdf(negModel.pdfs{i}, X_tmp)*negModel.inner_priors(i) ; 
        p1_z_c1 = horzcat(p1_z_c1, p_tmp) ; % p_x_giv_Cneg ;        
        p_tmp = p1_z_c0{i} + p1_z_c1{i} ;     
        p1_z = horzcat(p1_z, p_tmp) ;
        p_tmp = p1_z_c0{i} ./ (p1_z{i} + minTol) ;
        p1_c0_g_z = horzcat(p1_c0_g_z, p_tmp) ;
        
       
        sigmaPoints.X = horzcat(sigmaPoints.X, X_tmp) ;
        sigmaPoints.W = horzcat(sigmaPoints.W, W_tmp) ;  
    end
    precalcs.p1_z_c1 = p1_z_c1 ; % negative model evaluated
    precalcs.ppos = ppos ; % positive referece model evaluated posterior
    precalcs.p1_c0_g_z = p1_c0_g_z ;
    
    H.sigmaPoints = sigmaPoints ;
    H.precalcs = precalcs ;
    if isempty(posModel_r)
        return ;
    else
        negModel.precalcStat = H ;
    end
end
 
if isempty(negModel.precalcStat.sigmaPoints)
    H = 0 ;
    return ;
end

H_tmp = ones(1,length(negModel.precalcStat.precalcs.p1_z_c1)) ;

for i = 1 : length(negModel.precalcStat.precalcs.p1_z_c1)    
    
    X = negModel.precalcStat.sigmaPoints.X{i} ;
    W = negModel.precalcStat.sigmaPoints.W{i} ;
    p2_z_c0 = evaluatePointsUnderPdf(posModel_r, X)*modelPriors.pPos ; % p_x_Cpos
    
    p1_c0_g_z = negModel.precalcStat.precalcs.p1_c0_g_z{i} ;
    p1_c1_g_z = 1 - p1_c0_g_z ; % ; negModel.precalcStat.precalcs.p1_c1_g_z{i} ;  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
    p2_z_c1 = negModel.precalcStat.precalcs.p1_z_c1{i} ; %p1_z_c1 ;
    p2_z = p2_z_c0 + p2_z_c1 ;
    p2_c0_g_z = p2_z_c0 ./ (p2_z + minTol) ;
    p2_c1_g_z = 1 - p2_c0_g_z ; %p2_z_c1 ./ (p2_z + minTol) ;
 
    if type_cost == 1
        g0 = (sqrt(p1_c0_g_z) - sqrt(p2_c0_g_z)).^2 ;
        g1 = (sqrt(p1_c1_g_z) - sqrt(p2_c1_g_z)).^2 ;
        Hc = sum(W.*(g0 + g1)) ; H_tmp(i) = sqrt(Hc/2) ; 

%         Hc = sum(W.*sqrt((g0 + g1)/2)) ; H_tmp(i) =  Hc ; 
    else
        Hc = max(abs((p1_c0_g_z)-(p2_c0_g_z))) ; %         
        H_tmp(i) = Hc ;
    end
end
 
% figure(3); bar(sort(H_tmp,'descend')) ; 
% title(num2str(length(H_tmp))) ;
% axis([1,26,0,0.5]) ; drawnow ; pause(0.1) ; figure(1);

% H = sum(H_tmp) * (length(negModel.precalcStat.precalcs.p1_z_c1)+1) ; 
% H = sum(H_tmp)   * (length(negModel.precalcStat.precalcs.p1_z_c1)+1) ;
H = max(H_tmp) ; %* ((length(negModel.precalcStat.precalcs.p1_z_c1))) ;
 
% length(negModel.precalcStat.precalcs.p1_z_c1)
if isempty(H)
    H = 0 ;
end
 
% ----------------------- last good and working ------------------------ %
function H = uCostModel_prev( negModel, posModel, posModel_r, negModelPrior, use_approximate_calc, pdf_for_sigmas, type_cost )
% Matej Kristan (2009)
% calculates a cost of model reduction in terms of classification
% accuracy.
 

% od testa za Letter se ta razlikuje po tem, da sem aktiviral brisanje
% predhranjenja in doloèil v funkciji ki klièe to, da naj raèuna na toèkah
% s subseta, ki ga želimo kompresirati namesto na celotni pozitivni pdf.

% focus only on centers
% focus_on_centers = 0 ; % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

% use_approximate_calc = 1 ;
use_sigma_points = 0 ;

% if use_sigma_points == 1 
%     focus_on_centers = 0 ;
% else
%     focus_on_centers = 1 ;
% end

% check if negative model even exists for the local compression:
emptynegative = 0 ;
if isempty(negModel)  
    emptynegative = 1 ;
end
if emptynegative == 1 
   
    if isempty(posModel_r)
        H.sigmaPoints = [] ;
        H.precalcs = [] ;
    else
         H = 0 ;
    end
    return ;
end 

% negModel.precalcStat = [] ;

modelPriors.pNeg = negModelPrior ; % 0.5 ;
modelPriors.pPos = 1 - modelPriors.pNeg ; 

modelPriors.pPos = 0.5 ;
modelPriors.pNeg = 0.5 ;
% kasneje še negModel.inner_priors(i)=0.5 !!!!!!!!!!!!!!!!!!!!!!!!!

minTol = 1e-50 ;
nimPerc = 0.001 ;

for i = 1 : length(negModel.pdfs)
        negModel.inner_priors(i) = 0.5 ;
end

 
if ~isfield(negModel,'precalcStat') || isempty(negModel.precalcStat)
   % generate sigma points 
    MaxV = 3 ; 
    f0 = posModel ;
    
    % calculate sigma points for entire set
    if use_sigma_points == 1
       [X, sigPointsPerComponent, w, k ] = getAllSigmaPointsOnMixture( f0, MaxV ) ;        
        X = real(X) ;
        W = repmat(f0.w,sigPointsPerComponent,1) ;
        W = reshape(W,1,length(f0.w)*sigPointsPerComponent) ;
        w2 = repmat(w,1,length(f0.w)) ;
        W = W.*w2 ;
    else
        X = f0.Mu ;
        W = f0.w ;
    end
    % store vectors
    sigmaPoints.X = X ;
    sigmaPoints.W = W ;
    
 
    ppos = evaluatePointsUnderPdf(posModel, X) ;
       
    % evaluate probabilities of model 1
    p1_z_c0 = ppos*modelPriors.pPos ; % p_x_Cpos ;    
    p1_z_c1 = {} ; p1_z = {} ; p1_c0_g_z = {} ; %p1_c1_g_z = {} ;
    for i = 1 : length(negModel.pdfs)        
        negModel.inner_priors(i) = 0.5 ;        
        p_tmp = evaluatePointsUnderPdf(negModel.pdfs{i}, X)*negModel.inner_priors(i) ; 
        p1_z_c1 = horzcat(p1_z_c1, p_tmp) ; % p_x_giv_Cneg ;        
        p_tmp = p1_z_c0 + p1_z_c1{i} ;     
        p1_z = horzcat(p1_z, p_tmp) ;
        p_tmp = p1_z_c0 ./ (p1_z{i} + minTol) ;
        p1_c0_g_z = horzcat(p1_c0_g_z, p_tmp) ;
    end
    precalcs.p1_z_c1 = p1_z_c1 ; % negative model evaluated
    precalcs.p1_c0_g_z = p1_c0_g_z ; % positive referece model evaluated posterior

    H.sigmaPoints = sigmaPoints ;
    H.precalcs = precalcs ;
    if isempty(posModel_r)
        return ;
    else
        negModel.precalcStat = H ;
    end
end
 
if isempty(negModel.precalcStat.sigmaPoints)
    H = 0 ;
    return ;
end

H_tmp = ones(1,length(negModel.precalcStat.precalcs.p1_z_c1)) ;
X = negModel.precalcStat.sigmaPoints.X ;
W = negModel.precalcStat.sigmaPoints.W ;
p2_z_c0 = evaluatePointsUnderPdf(posModel_r, X)*modelPriors.pPos ; % p_x_Cpos
for i = 1 : length(negModel.precalcStat.precalcs.p1_z_c1)    
    p1_c0_g_z = negModel.precalcStat.precalcs.p1_c0_g_z{i} ;
    p1_c1_g_z = 1 - p1_c0_g_z ; % ; negModel.precalcStat.precalcs.p1_c1_g_z{i} ;  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
    p2_z_c1 = negModel.precalcStat.precalcs.p1_z_c1{i} ; %p1_z_c1 ;
    p2_z = p2_z_c0 + p2_z_c1 ;
    p2_c0_g_z = p2_z_c0 ./ (p2_z + minTol) ;
    p2_c1_g_z = 1 - p2_c0_g_z ; %p2_z_c1 ./ (p2_z + minTol) ;
 
    
    if type_cost == 1
        g0 = ((sqrt(p1_c0_g_z) - sqrt(p2_c0_g_z)).^2) ;
        g1 = ((sqrt(p1_c1_g_z) - sqrt(p2_c1_g_z)).^2) ;
        Hc = sum(W.*(g0 + g1)) ;
        H_tmp(i) = sqrt(Hc/2) ;          
    else
        Hc = max(abs((p1_c0_g_z)-(p2_c0_g_z))) ; %         
        H_tmp(i) = Hc ;
    end
end

H = max(H_tmp) ;
if isempty(H)
    H = 0 ;
end

