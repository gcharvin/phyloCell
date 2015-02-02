function H = localDiscCostModel( negModel, posModel, posModel_r, negModelPrior, use_approximate_calc )
% Matej Kristan (2009)
% calculates a cost of model reduction in terms of classification
% accuracy.
 
% use_approximate_calc = 1 ;

modelPriors.pNeg = negModelPrior ;
modelPriors.pPos = 1 - modelPriors.pNeg ;
% modelPriors.pNeg = 0.5 ;
% modelPriors.pPos = 1 - modelPriors.pNeg ;
minTol = 1e-50 ;
nimPerc = 0.01 ;
 
% pdf1.z_g_c0 = posModel ;
% pdf1.z_g_c1 = negModel ;
% pdf1.c01 = [modelPriors.pPos, modelPriors.pNeg ] ;
% 
% pdf2.z_g_c0 = posModel_r ;
% pdf2.z_g_c1 = negModel ;
% pdf2.c01 = [modelPriors.pPos, modelPriors.pNeg ] ;
 
if 1 ==1   || ~isfield(negModel,'precalcStat') || isempty(negModel.precalcStat)
   % generate sigma points 
    MaxV = 3 ;
    f0 = mergeDistributions( posModel, posModel_r, [0.5, 0.5] ) ;
%     f0 = posModel ;
    
    % calculate sigma points for entire set
    [X, sigPointsPerComponent, w, k ] = getAllSigmaPointsOnMixture( f0, MaxV ) ;
    X = real(X) ;
    W = repmat(f0.w,sigPointsPerComponent,1) ;
    W = reshape(W,1,length(f0.w)*sigPointsPerComponent) ;
    w2 = repmat(w,1,length(f0.w)) ;
    W = W.*w2 ;
 
    % for approximate calculation
    if use_approximate_calc == 1
         % determine the relevant region 
         id_lastX = length(posModel.w)*sigPointsPerComponent ;
         minX = min(X(:,1:id_lastX), [], 2) ;
         maxX = max(X(:,1:id_lastX), [], 2) ;
         new_mu = (minX + maxX)/2 ;
         new_Cov = diag((new_mu - minX).^2 + 1e-20)  ;
         pmax = normpdf(new_mu, new_mu, [], new_Cov) ;
         peval = normpdf(X, new_mu, [], new_Cov) ;
         id_sel = peval / pmax > nimPerc ; % 0.01 ; %0.2 ;
%          tt = sum(id_sel) ;
%         
%         [new_mu, new_Cov, w_out] = momentMatchPdf(posModel.Mu, posModel.Cov, posModel.w) ;
%         pmax = normpdf(new_mu, new_mu, [], new_Cov) ;
%         peval = normpdf(X, new_mu, [], new_Cov) ;
%         id_sel = peval / pmax > 1e-3 ;
%         sum(id_sel) - tt
        
        X = X(:,id_sel) ;
        W = W(:,id_sel) ;
        W = W / sum(W) ;
    end
 
    % store vectors
    sigmaPoints.X = X ;
    sigmaPoints.W = W ;
    
    ppos = evaluatePointsUnderPdf(posModel, X) ;
%     f0_val = evaluatePointsUnderPdf(f0, X) ;  
%     sigmaPoints.W = W.*ppos ./ f0_val ;
    
    
    % evaluate probabilities of model 1
    p1_z_c0 = ppos*modelPriors.pPos ;
    p1_z_c1 = evaluatePointsUnderPdf(negModel, X)*modelPriors.pNeg ;
    p1_z = p1_z_c0 + p1_z_c1 ; 
    p1_c0_g_z = p1_z_c0 ./ (p1_z + minTol) ;
    p1_c1_g_z = p1_z_c1 ./ (p1_z + minTol) ;
 
%     % for approximate calculation
%     if use_approximate_calc == 1
% %         id_sel = p1_c0_g_z > 1e-5 ;         
%         sigmaPoints.X = X(:,id_sel) ;
%         sigmaPoints.W = W(:,id_sel) ;
%         precalcs.p1_z_c1 = p1_z_c1(:,id_sel)  ;
%         precalcs.p1_c0_g_z = p1_c0_g_z(:,id_sel)  ;
%         precalcs.p1_c1_g_z = p1_c1_g_z(:,id_sel)  ;
%     else        
        precalcs.p1_z_c1 = p1_z_c1 ;
        precalcs.p1_c0_g_z = p1_c0_g_z ;
        precalcs.p1_c1_g_z = p1_c1_g_z ;
%     end
    
    H.sigmaPoints = sigmaPoints ;
    H.precalcs = precalcs ;
    if isempty(posModel_r)
        return ;
    else
        negModel.precalcStat = H ;
    end
end
 

p1_z_c1 = negModel.precalcStat.precalcs.p1_z_c1 ;
p1_c0_g_z = negModel.precalcStat.precalcs.p1_c0_g_z ;
p1_c1_g_z = negModel.precalcStat.precalcs.p1_c1_g_z ;
X = negModel.precalcStat.sigmaPoints.X ;
W = negModel.precalcStat.sigmaPoints.W ;
 
% calculate error value
p2_z_c1 = p1_z_c1 ;
p2_z_c0 = evaluatePointsUnderPdf(posModel_r, X)*modelPriors.pPos ;
p2_z = p2_z_c0 + p2_z_c1 ; %p1_z ;%5
p2_c0_g_z = p2_z_c0 ./ (p2_z + minTol) ;
p2_c1_g_z = p2_z_c1 ./ (p2_z + minTol) ;

% classification error
% P2_clsf1 = p2_c0_g_z >= 0.5 ;
% P1_clsf1 = p1_c0_g_z >= 0.5 ;
% D_clsf = abs(P2_clsf1 - P1_clsf1) ;
% H = sum(W.*D_clsf) ;


% % % hellinger over posterior
g0 = ((sqrt(p1_c0_g_z) - sqrt(p2_c0_g_z)).^2) ;
g1 = ((sqrt(p1_c1_g_z) - sqrt(p2_c1_g_z)).^2) ;
Hc = sum(W.*(g0 + g1)) ;
H = sqrt(Hc/2) ;

% g0 = abs(((p1_c0_g_z) - (p2_c0_g_z))) ;
% g1 = abs(((p1_c1_g_z) - (p2_c1_g_z))) ;
% H = max([g0,g1]) ; 
 

% H = sum(W.*sqrt((g0 + g1)/2)) ;