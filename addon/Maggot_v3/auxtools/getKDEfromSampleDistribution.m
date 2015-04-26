function model = getKDEfromSampleDistribution( model, N_eff ) 
 
% model.smod.useVbw = 0 ;
N = size(model.w,2) ;
% model.Cov = {} ;
if isempty(model.Cov)
    for i = 1 : N
        model.Cov = horzcat(model.Cov, model.smod.ps.Cov{i}+model.smod.H) ;
    end
else
    for i = 1 : N 
        model.Cov{i} = model.smod.ps.Cov{i} + model.smod.H ;
    end
    
% %     if isfield(model.smod,'useVbw') && model.smod.useVbw == 1
%         a = - 2/(2*(size(model.Mu,1) + 4)) ;
%         for i = 1 : N
%             alph = ((N_eff*model.w(i))^a) ;
%             w0 = 1 + alph*(1-alph) ; w1 = alph ;
%             model.Cov{i} = model.smod.ps.Cov{i}*w0 + model.smod.H*w1 ;
%         end
% %     else
% %     for i = 1 : N
% %         model.Cov{i} = model.smod.ps.Cov{i} + model.smod.H ;
% %     end
end
 