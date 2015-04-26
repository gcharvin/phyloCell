%%
% Originally a part of: Maggot (developed within EU project CogX)
% Author: Matej Kristan, 2009 (matej.kristan@fri.uni-lj.si; http://vicos.fri.uni-lj.si/matejk/)
% Last revised: 2009
%%
function  [H_opt, F_opt, h_amise, dim_subspace] = leastSquaresCrossValidation( varargin )
% global pdf_lscv ;

minValAllowed = 1e-5 ; % was 1e-10
use_removeNullspace = 0 ;
likth = 1e-5 ;
% ikdeParams = [] ;
minTol = 0.01 ; % 1e-3 ;
maxNumIt = 50 ; % 20
model = [] ;
kernelType = 'general' ; % 'spherical'
obs = [] ;
% process arguments
args = varargin;
nargs = length(args);
for i = 1:2:nargs
    switch args{i}   
        case 'model', model = args{i+1} ;
        case 'obs', obs = args{i+1} ;
        case 'kernelType', kernelType = args{i+1} ; 
        case 'ikdeParams', ikdeParams = args{i+1} ; 
    end
end

if ~isempty(model)
    if ~isempty(model.Mu)
        error('LSCV does not support using mixture models for input data yet!') ;
%         kernelType = 'spherical' ; %         warning('Switching to spherical kernel.') ;
    end
end
 
% if use_removeNullspace == 1
%     derivativeModel = augmentMixtureModelByThis( model, obs,...
%         0, ikdeParams.obs_mixing_weights,...
%                         ikdeParams.mix_weights ) ;  
%     [derivativeModel, H_o ]= readjustKernels( derivativeModel, 0 ) ;
% 
%     minEigenEnergy = 1e-5 ;
%     output = subspacePrewhitenTransform( 'pdf', derivativeModel, 'globalCov',...
%         ikdeParams.scale.Cov, 'minEigenEnergy', minEigenEnergy, ...
%         'transDirection', 'forward',...
%         'allLayers', 0 ) ;
%     pdft = output.pdf ;
%     obs = pdft.Mu ;
%     globalCov = output.globalCov ;% eye(size(output.globalCov)) ;%output.globalCov ;
%     % warning('Prewhiten using moments!! not sample cov!!')
%     invF_trns = output.svdRes.V*sqrt(abs(output.svdRes.S)) ;
%     d = size(pdft.Mu,1) ;
%     dim_subspace = d ;
%     C0 = globalCov ;
% else
    C0 = cov(obs') ;
% end

[U,S,V] = svd(C0) ;
% s = diag(S) ; s(s<minValAllowed) = min(s(s>=minValAllowed)) ;
% S = diag(s) ;
% 
 min_s = 1e-4 ; 
    s = diag(S) ; ss = s/max(s) ; s(ss<min_s) = min(s(ss>min_s)); 
    S = diag(s) ; minValAllowed = 0;   


C0 = U*S*U' ;

N = size(obs,2) ;
% Ct_1 = C0 *((4/(3*N))^(1/5))^2 * 0.8^2 ; %0.85^2 ;
 d = size(obs,1);  Ct_1 = C0 *(4/((d+2)*N))^(2/(d+4)) * 0.8^2 ;
dim_subspace = size(Ct_1,1)  ;
min_err = 1e-30 ;
loglik_prev = [] ;
switch kernelType
    case 'general'        
        Ct = Ct_1*0 ;
        idx = 1 : N ;
        iter_cnt = 0 ; 
        while iter_cnt < maxNumIt
            I = 0 ;  
            Num_sumterms = 0 ;
            loglik = 1 ;
            for i = 1 : N
                s_idx = idx ~= i ;
%                 Gij = normpdf(obs(:,s_idx), obs(:,i), [], Ct_1) ;
                Gij = normpdf(obs(:,s_idx), obs(:,i), [], Ct_1) ;
 
                if isnan(sum(sum(Gij)) )
                    sdg
                end
                
                S = 0 ;
                i_cnt = 1 ;
                for j = 1 : N
                    if (i==j)
                        continue ;
                    end                    
                    
                    d = obs(:,i) -  obs(:,j) ;
                    S = S + d*d'* Gij(i_cnt) ; %normpdf(obs(:,j), obs(:,i), [], Ct_1) ;%
                    i_cnt = i_cnt + 1 ;
                     Num_sumterms = Num_sumterms + 1 ; 
                end
                if isnan(sum(sum(S)) )
                    sdg
                end
%                 I = I + S/(sum(Gij)/(N-1)) ;
                I = I + S / max([mean(Gij),min_err]) ;
                if isnan(sum(sum(I)) )
                    sdg
                end
                loglik = loglik + -log(max([mean(Gij),min_err])) ;
            end
            Ct = I / Num_sumterms ; %( N*(N-1) ) ;
            loglik = loglik / Num_sumterms ;
            [Uct,Sct,Vct] = svd(Ct) ;
            
              s = diag(Sct) ; ss = s/max(s) ; s(ss<min_s) = min(s(ss>min_s)); 
              Sct = diag(s) ; 
            
            dSct = diag(Sct) ;    
            id_critical = sum(dSct < minValAllowed) ;
            dSct(dSct < minValAllowed) = min(dSct(dSct >= minValAllowed)) ;  
            Sct = diag(dSct) ;
            
 
            
            Ct = Uct*Sct*Uct' ;
%             v = abs(Ct - Ct_1) ; v = v(:)' ;
%             vd = abs(Ct_1) ; vd = max([vd(:)'*0+1e-32; vd(:)']) ;
%             v = v./vd ;
%             if mean(v) < minTol
%                 break ;
%             end            
%             vd = max([diag(Ct_1)';diag(Ct_1)'*0 + 1e-32])' ;
%             delta = mean(abs(diag(Ct)-diag(Ct_1))./vd) ;
%             if delta < minTol
%                 break ;
%             end
            if iter_cnt == 0 
                loglik_prev = loglik ;
            else
                if id_critical == 1 || loglik < 0 || loglik > loglik_prev
                    Ct = Ct_1 ;
                    break ;
                end
                
                if abs(loglik - loglik_prev) < likth %abs(loglik - loglik_prev)/(loglik_prev) < minTol
                    break ;
                end
            end

            Ct_1 = Ct ;
            iter_cnt = iter_cnt + 1 ;
            loglik_prev = loglik ;
        end
        if iter_cnt == maxNumIt
            warning('LSCV: Maximum number of iteration reached!') ;
        end
    case 'spherical'
%         pdf_lscv = pdft ;
%         c0 = 1e-20 ;
%         [ F_c1, c1 ] = getStructFromBW( Ct_1/0.7^2 ) ;
%         c = fminbnd(@(x)costLogLik(x),c0, c1) ;
%         Ct = eye(size(pdft.Mu,1))*c ;
    error('Not implemented yet!') ;
end

H_opt = Ct  ;
% [ F_opt, h_amise ] = getStructFromBW( H_opt ) ;

if use_removeNullspace == 1
    % augment H by nullspace
    C_prot = zeros(size(output.svdRes.S)) ;
    C_prot(output.svdRes.id_valid,output.svdRes.id_valid) = H_opt ;
    H_opt = C_prot ;

    % backtransform bandwidth and calculate the structure
    H_opt = invF_trns*H_opt*invF_trns' ;
end
[ F_opt, h_amise ] = getStructFromBW( H_opt ) ;



function gL = costLogLik( c )
global pdf_lscv ;
% C1 = C0 + c*I ;

N = length(pdf_lscv.w) ;
d = size(pdf_lscv.Mu,1) ;
I = eye(d) ;
C = c*I ;
gL = 0 ;
for i = 1 : N
    p_mui = 0 ;
    gL_inner = 0 ;
    iNi = 0 ;
    for j = 1 : N
        if (i==j)
            continue ;
        end
        Cij = pdf_lscv.Cov{i} + pdf_lscv.Cov{j} + C ;
        invCij = inv(Cij) ;
        D = pdf_lscv.Mu(:,i) - pdf_lscv.Mu(:,j) ;
        Gij = normpdf(D,D*0,[],Cij) ;
        gL_inner = gL_inner + pdf_lscv.w(j)*((invCij*D*D' - I)*invCij - invCij*Gij) ;
        iNi = iNi + pdf_lscv.w(j) ;
        p_mui = pdf_lscv.w(j)*Gij ;
    end
    Ni = 1 / iNi ;
    p_mui = p_mui * Ni ;
    gL = gL + pdf_lscv.w(i)*(1/p_mui)*gL_inner ;    
end
gL = det(gL) ;



     