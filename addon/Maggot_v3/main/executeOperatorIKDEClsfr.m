function hyper_output_kde_cl = executeOperatorIKDEClsfr( hyper_input_kde_cl, varargin )
% If inputs are data-points
% input_data{i}.class ... with data we need to specify also the class label
% input_data{i}.data  ... actual input data-points

% If inputs are kdes
% input_data{i}.class ... class label
% input_data{i}.kde   ... kde

% by default the classifier construction is discriminative rather than reconstructive


%  hyper_input_kde_cl.kde_cl = classes_kde_cl ;
%  hyper_input_kde_cl.typeRecDescr = 'discriminative' ; % {discriminative, reconstructive}
%  hyper_input_kde_cl.class_labels = class_labels ;
%

Cdescribe = 3 ;
Cans = 1 ; 0.1*0 ; 
Clist = 0.1*0 ;

turn_off_splitting = [] ;
type_update = 'partial' ; % 'partial', 'joint'
% useSomeOtherKindOfEstimator = [] ; % empty means default, 'adaptiveMixtures' uses adaptive mixtures
use_equalimportance =  0 ;
classifyWithPosterior = 1 ;
extensive_answer = 0 ;
minNumDataPointsToFormKDE = [] ; 
autonomous_update = [] ; 
exclude_model = [] ;

force_value_init_of_maxNumCompsBeforeCompression = [] ;
features_to_show = [] ;
pair_dist_struct_use_approx = [] ;
random_fselect_threshold = [] ;
compressionClusterThresh = [] ;
autoUpdateThres_upper = [] ;
autoUpdateThres_lower = [] ;
autonomous_update = [] ; %  'oracle_verified', 'self_verified' updates only if it misclassifies the input data, 'pure_oracle'
% not operational % minimal_required_examps_mode = [] ; % if 1, then a vector is added only if it's wrongly classified
min_samps_per_model_feat_sel = [] ;
react_compression_to_feature_selection = [] ;
min_th_feat_sel = [] ;
clear_selected_features = 0 ;
sub_selected_features = [] ;
val_set = [] ;
val_get = [] ;
use_unknown_model = [] ;
unknown_model_value = [] ;
create_clsfr_from_this = [] ;
ignoreClasses = [] ;
typeRecDescr = [] ; % type of model: {discriminative, reconstructive}
operator_data = ' ' ;
input_data = [] ;
vforwvargin = {} ;
args = varargin;
nargs = length(args);
i = 1 ;
while i <= nargs   
    switch args{i}                
        case 'input_data', input_data = args{i+1} ; i = i + 2 ; 
        case 'make_simple_feature_selection', operator_data = args{i} ; i = i + 1 ;   
        case 'init', operator_data = args{i} ; i = i + 1 ;
        case 'auto_choose_model', operator_data = args{i} ; i = i + 1 ;     
        case 'autonomous_update', autonomous_update = args{i+1} ; i = i + 2 ;
        case 'get_class_names', operator_data = args{i} ; i = i + 1 ;
        case 'get_class_indexes', operator_data = args{i} ; i = i + 1 ;  
        case 'get_name_at_index', operator_data = args{i} ; val_get = args{i+1} ; i = i + 2 ;            
        case 'get_index_at_name', operator_data = args{i} ; val_get = args{i+1} ; i = i + 2 ; 
        case 'set_name_at_index', operator_data = args{i} ; val_get = args{i+1} ; val_set = args{i+2} ; i = i + 3 ;         
        case 'add_input', operator_data = args{i} ; i = i + 1 ; 
        case 'introspect', operator_data = args{i} ; i = i + 1 ;     
        case 'kdes_to_classifier', operator_data = args{i} ; i = i + 1 ;    
        case 'compress_pdf', operator_data = args{i} ; i = i + 1 ;
        case 'classifyData', operator_data = args{i} ; i = i + 1 ; 
        case 'classifyDataInSubspace', operator_data = args{i} ; i = i + 1 ;      
        case 'predict_missing', operator_data = args{i} ; i = i + 1 ; 
        case 'calculate_gains', operator_data = args{i} ; i = i + 1 ;
        case 'unlearn_with_input', operator_data = args{i} ; i = i + 1 ;
        case 'showKDE_of_class_index', operator_data = 'showKDE_of_class_index' ; val_get = args{i+1} ; i = i + 2 ;
        case 'regularize_all_kdes', operator_data = args{i} ; i = i + 1 ;
        case 'pair_dist_struct_use_approx', pair_dist_struct_use_approx = args{i+1} ; i = i + 2 ;
        case 'exclude_model', exclude_model = args{i+1} ; i = i + 2 ;
        case 'random_fselect_threshold', random_fselect_threshold = args{i+1} ; i = i + 2 ;
        case 'autoUpdateThres_upper', autoUpdateThres_upper = args{i+1} ; i = i + 2 ;    
        case 'autoUpdateThres_lower', autoUpdateThres_lower = args{i+1} ; i = i + 2 ;    
        case 'min_samps_per_model_feat_sel', min_samps_per_model_feat_sel = args{i+1} ; i = i + 2 ;
        case 'sub_selected_features', sub_selected_features = args{i+1} ; i = i + 2 ;
        case 'clear_selected_features', clear_selected_features = 1 ; i = i + 1 ;
        case 'min_th_feat_sel', min_th_feat_sel = args{i+1} ; i = i + 2 ;
        case 'react_compression_to_feature_selection', react_compression_to_feature_selection = args{i+1} ; i = i + 2 ;
        case 'typeRecDescr', typeRecDescr = args{i+1} ; i = i + 2 ;   
        case 'ignoreClasses', ignoreClasses = args{i+1} ; i = i + 2 ; 
        case 'classifyWithPosterior', classifyWithPosterior = args{i+1} ; i = i + 2 ; 
        case 'use_equalimportance', use_equalimportance = args{i+1} ; i = i + 2 ; 
        case 'type_update', type_update = args{i+1} ; i = i + 2 ;     
        case 'unknown_model_value', unknown_model_value = args{i+1} ; i = i + 2 ; 
        case 'use_unknown_model', use_unknown_model = args{i+1} ; i = i + 2 ; 
        case 'minNumDataPointsToFormKDE', minNumDataPointsToFormKDE = args{i+1} ; i = i + 2 ; 
        case 'compressionClusterThresh', compressionClusterThresh = args{i+1} ; i = i + 2 ;       
        case 'force_value_init_of_maxNumCompsBeforeCompression', force_value_init_of_maxNumCompsBeforeCompression = args{i+1} ; i = i + 2 ;       
%         case 'minimal_required_examps_mode', minimal_required_examps_mode
%         = args{i+1} ; i = i + 2 ;  % not operational
        case 'features_to_show', features_to_show = args{i+1} ; i = i + 2 ;           
        case 'extensive_answer', extensive_answer = 1 ; i = i + 1 ;         
        case 'useSomeOtherKindOfEstimator', 
            if ~isequal(args{i+1},0)
                vforwvargin = horzcat(vforwvargin, args{i} ) ;
                vforwvargin = horzcat(vforwvargin, args{i+1} ) ;                 
            end
            i = i + 2 ; 
        case 'compressionClusterThresh', 
            vforwvargin = horzcat(vforwvargin, args{i} ) ;
            vforwvargin = horzcat(vforwvargin, args{i+1} ) ;
            i = i + 2 ; 
        case 'showkdecolor'
             vforwvargin = horzcat(vforwvargin, args{i} ) ;
            vforwvargin = horzcat(vforwvargin, args{i+1} ) ;
            i = i + 2 ; 
        case 'turn_off_splitting', turn_off_splitting = args{i+1} ; i = i + 2 ;
        otherwise
             vforwvargin = horzcat(vforwvargin, args{i} ) ;
             i = i + 1 ; 
    end
end

% initialize parameters if required
if isempty(hyper_input_kde_cl)
    hyper_input_kde_cl.kde_cl = {} ;
    hyper_input_kde_cl.typeRecDescr = 'discriminative' ; % {discriminative, reconstructive}
    hyper_input_kde_cl.class_labels = [] ;
    hyper_input_kde_cl.class_labels_names = {} ;
    hyper_input_kde_cl.class_gains = [] ;
    hyper_input_kde_cl.unknown_model_value = 0 ;
    hyper_input_kde_cl.use_unknown_model = 0 ;
    hyper_input_kde_cl.Params.minNumDataPointsToFormKDE = 2 ;
    hyper_input_kde_cl.sub_selected_features = [] ;
    hyper_input_kde_cl.min_th_feat_sel = 0.1 ;
    hyper_input_kde_cl.min_samps_per_model_feat_sel = 5 ;
    hyper_input_kde_cl.react_compression_to_feature_selection = 0 ;
    hyper_input_kde_cl.sub_feature_sel_forgetting = 0 ;
    hyper_input_kde_cl.cummulative_feat_costs = [] ;
    hyper_input_kde_cl.autoUpdateThres.upper = 0.15 ;
    hyper_input_kde_cl.autoUpdateThres.lower = -1 ;    
    
    hyper_input_kde_cl.compressionClusterThresh.thReconstructive = 0.0500 ;
    hyper_input_kde_cl.compressionClusterThresh.thDiscriminative = 0.0500 ;
    hyper_input_kde_cl.random_fselect_threshold = 1 ;
    
    pair_dist_struct.dist = [] ;
    pair_dist_struct.dist_th = 4^2 ; %2.0^2 ;
    pair_dist_struct.use_approx = 1 ;
    hyper_input_kde_cl.pair_dist_struct = pair_dist_struct ;
    hyper_input_kde_cl.force_value_init_of_maxNumCompsBeforeCompression = -1 ;
    hyper_input_kde_cl.turn_off_splitting = 0 ;
%     hyper_input_kde_cl.minimal_required_examps_mode = 0 ;
end

if ~isempty(turn_off_splitting)
    hyper_input_kde_cl.turn_off_splitting = turn_off_splitting ;
end
 
if ~isempty(force_value_init_of_maxNumCompsBeforeCompression)
    hyper_input_kde_cl.force_value_init_of_maxNumCompsBeforeCompression = force_value_init_of_maxNumCompsBeforeCompression ;
end

if ~isempty(pair_dist_struct_use_approx)
    hyper_input_kde_cl.pair_dist_struct.use_approx  = pair_dist_struct_use_approx ;
end
 
if ~isempty(random_fselect_threshold)
    hyper_input_kde_cl.random_fselect_threshold  = random_fselect_threshold ;
end

if ~isempty(compressionClusterThresh)
    hyper_input_kde_cl.compressionClusterThresh  = compressionClusterThresh ;
end
 
if ~isempty(min_samps_per_model_feat_sel)
    hyper_input_kde_cl.min_samps_per_model_feat_sel  = min_samps_per_model_feat_sel ;
end

if ~isempty(react_compression_to_feature_selection)
    hyper_input_kde_cl.react_compression_to_feature_selection  = react_compression_to_feature_selection ;
end

if ~isempty(min_th_feat_sel)
    hyper_input_kde_cl.min_th_feat_sel  = min_th_feat_sel ;
end
 
if clear_selected_features ~= 0
    hyper_input_kde_cl.sub_selected_features = [] ;
end

if ~isempty(sub_selected_features)
   hyper_input_kde_cl.sub_selected_features = sub_selected_features ; 
end

if ~isempty(minNumDataPointsToFormKDE)
    hyper_input_kde_cl.Params.minNumDataPointsToFormKDE = minNumDataPointsToFormKDE ;
end

if ~isempty(typeRecDescr)
    hyper_input_kde_cl.typeRecDescr = typeRecDescr ;
    
    if isequal(typeRecDescr, 'dKDE')
          hyper_input_kde_cl.turn_off_splitting = 1 ; 
    end
    
end

if ~isempty(unknown_model_value)
    hyper_input_kde_cl.unknown_model_value = unknown_model_value ;
end

if ~isempty(use_unknown_model)
   hyper_input_kde_cl.use_unknown_model = use_unknown_model  ;
end

if ~isempty(autonomous_update)           
    if isequal(autonomous_update, 'pure_oracle')           
       operator_data = 'add_input' ;      
    else
        operator_data = 'autonomous_update' ; 
    end
end

if ~isempty(autoUpdateThres_upper)
    hyper_input_kde_cl.autoUpdateThres.upper = autoUpdateThres_upper ;
end

if ~isempty(autoUpdateThres_lower)
    hyper_input_kde_cl.autoUpdateThres.lower = autoUpdateThres_lower ;
end
 
switch operator_data    
    case 'init'
        hyper_output_kde_cl = hyper_input_kde_cl ;
    case 'get_name_at_index'
        hyper_output_kde_cl = [] ;
        if isempty(val_get)
            return ;
        end
        
        if val_get == -1 || val_get == 0
            hyper_output_kde_cl = '0' ; % 'unknown'
            return ;
        end
        
        if val_get <= length(hyper_input_kde_cl.class_labels)
            hyper_output_kde_cl = hyper_input_kde_cl.class_labels_names{val_get} ;
        else
            return ;
        end
    case 'get_index_at_name' 
        hyper_output_kde_cl = [] ;
        if isempty(val_get)
            return ;
        end
        
        for i = 1 : length(hyper_input_kde_cl.class_labels_names)
           if isequal(val_get, hyper_input_kde_cl.class_labels_names{i}) 
               hyper_output_kde_cl = i ;
               return ;
           end
        end
        return ;
    case 'set_name_at_index'
        hyper_output_kde_cl = hyper_input_kde_cl ;
        if isempty(val_set) || val_get > length(hyper_output_kde_cl.class_labels)
            return ;
        end
        hyper_output_kde_cl.class_labels_names{val_get} = val_set;
        return ;        
    case 'get_class_indexes'
        hyper_output_kde_cl = hyper_input_kde_cl.class_labels ;
        return ;
    case 'get_class_names'
        hyper_output_kde_cl = hyper_input_kde_cl.class_labels_names ;
        return ;    
    case 'kdes_to_classifier'   
        class_labels = [] ;
        classes_kde_cl = {} ;
        class_labels_names = {} ;
        for i = 1 : length(input_data)
            
            % initialize cumulative feature selector
            if isempty(hyper_input_kde_cl.cummulative_feat_costs)
               hyper_input_kde_cl.cummulative_feat_costs = zeros(1,size(input_data{i}.kde.pdf.Mu,1)) ;  
            end
            
            classes_kde_cl = horzcat(classes_kde_cl, input_data{i}.kde) ;
            class_labels = horzcat(class_labels, input_data{i}.class) ;
            class_labels_names = horzcat(class_labels_names, num2str(input_data{i}.class)) ;
        end        
        hyper_input_kde_cl.kde_cl = classes_kde_cl ;
        hyper_input_kde_cl.typeRecDescr = 'discriminative' ; % {discriminative, reconstructive}
        hyper_input_kde_cl.class_labels = class_labels ;
        
        % recalculate the unknown model value!!
        hyper_input_kde_cl.unknown_model_value = prob_of_unknown( hyper_input_kde_cl.kde_cl ) ;
        hyper_output_kde_cl = hyper_input_kde_cl ;
    case 'make_simple_feature_selection'
        r = rand(1) ; 
        if r > hyper_input_kde_cl.random_fselect_threshold
            hyper_output_kde_cl = hyper_input_kde_cl ;
            return ;
        end
        
        % check if there is enough samples per model
        N_eff = 0 ;
        for i = 1 : length(hyper_input_kde_cl.kde_cl)
             N_eff = N_eff + hyper_input_kde_cl.kde_cl{i}.ikdeParams.N_eff ;             
        end
        if  N_eff < hyper_input_kde_cl.min_samps_per_model_feat_sel
            hyper_output_kde_cl = hyper_input_kde_cl ;
            return ;
        end
 
        d_before = hyper_input_kde_cl.sub_selected_features ;
        fast_search = 1 ;
        % prepare distributions  
        pdf_classes.N_data = 0 ;
        pdf_classes.F_giv_cj = {} ;
        pdf_classes.cj = [] ;
        for i = 1 : length(hyper_input_kde_cl.kde_cl)
%             [new_mu, new_Cov, w_out] = momentMatchPdf(hyper_input_kde_cl.kde_cl{i}.pdf.Mu, hyper_input_kde_cl.kde_cl{i}.pdf.Cov, hyper_input_kde_cl.kde_cl{i}.pdf.w) ;
%             p.Mu = new_mu ;
%             p.Cov = {new_Cov} ;
%             p.w = 1 ;
           % get kde from sublayer
%             pdf_tmp = extractMixtureFromSublayers( hyper_input_kde_cl.kde_cl{i}.pdf ) ;

            p.Mu = hyper_input_kde_cl.kde_cl{i}.pdf.Mu ; %pdf_tmp.Mu ; % 
            p.Cov = hyper_input_kde_cl.kde_cl{i}.pdf.Cov ; %pdf_tmp.Cov ; % 
% %             for i_cv = 1 : length(p.Cov)
% %                p.Cov{i_cv} = regularizeCovariance( p.Cov{i_cv} )  ;
% %                p.Cov{i_cv} = p.Cov{i_cv} + eye(size(p.Cov{i_cv}))*(1e-5) ;
% %             end
            p.w = hyper_input_kde_cl.kde_cl{i}.pdf.w ; %pdf_tmp.w ; %
            pdf_classes.F_giv_cj = horzcat(pdf_classes.F_giv_cj, p) ;
            pdf_classes.N_data = pdf_classes.N_data + hyper_input_kde_cl.kde_cl{i}.ikdeParams.N_eff ;
            pdf_classes.cj = [pdf_classes.cj, hyper_input_kde_cl.kde_cl{i}.ikdeParams.N_eff] ;
        end
        pdf_classes.cj = pdf_classes.cj / sum(pdf_classes.cj) ;
  
        
          % select features                
        [ Cost, f_sel ] = mixtureFeatureSelection( pdf_classes, fast_search ) ;  
        id_cost = 2 ;
        rslt = Cost(3:size(Cost,1),:) ;
%         
%         max_v = (rslt(id_cost,size(rslt,2))) ;
%         
%         

        indic = rslt(id_cost,:) > hyper_input_kde_cl.min_th_feat_sel*0.1 ; %max_v*hyper_input_kde_cl.min_th_feat_sel ;
        frst = 1 ; 
        for i = length(indic) :-1: 1
           if indic(i) == 0 
               frst = i + 1;
               break ;
           end
        end
        indic = indic*0 ;
        indic(frst:length(indic)) = 1 ;
        sel_feats = sort(rslt(1, find(indic))) ;
        hyper_input_kde_cl.sub_selected_features = sel_feats ;      
        hyper_output_kde_cl = hyper_input_kde_cl ;
        
        hyper_output_kde_cl.debug.Cost = Cost ;
        
        [srtd, resuf_idxs]=sort(Cost(3,:)) ;
        
        if sum(hyper_output_kde_cl.cummulative_feat_costs) == 0 
            hyper_output_kde_cl.cummulative_feat_costs = rslt(3,resuf_idxs)/max(rslt(3,resuf_idxs)) ;
        else
            weigh = 0.9 ;
            rr = rslt(3,resuf_idxs)/max(rslt(3,resuf_idxs)) ;
            hyper_output_kde_cl.cummulative_feat_costs =  hyper_output_kde_cl.cummulative_feat_costs*weigh + (1-weigh) * rr ;
        end
%         ratio = hyper_output_kde_cl.cummulative_feat_costs/max(hyper_output_kde_cl.cummulative_feat_costs)
%         hyper_input_kde_cl.sub_selected_features = find( ratio > 0.2 ) ;
        
%         hyper_input_kde_cl.sub_selected_features = find(hyper_output_kde_cl.cummulative_feat_costs > hyper_input_kde_cl.min_th_feat_sel ) ; 
%         rslt(id_cost,rslt(1,:))
%         hyper_input_kde_cl.sub_selected_features
  hyper_input_kde_cl.sub_selected_features = find(hyper_output_kde_cl.cummulative_feat_costs > 0.65) ;
  [rslt(2,resuf_idxs);rslt(3,resuf_idxs)/max(rslt(3,resuf_idxs));hyper_output_kde_cl.cummulative_feat_costs ]
  hyper_input_kde_cl.sub_selected_features

%         [ix, xsrt] = sort(rslt(1,:)) ;
%         cumad = rslt(id_cost,:).*(indic+ 0.05*(indic==0)) ;
%         hyper_output_kde_cl.cummulative_feat_costs = hyper_output_kde_cl.cummulative_feat_costs + cumad(xsrt) ; 
%         hyper_output_kde_cl.cummulative_feat_costs
% %         % select features                
% %         [ Cost, f_sel ] = mixtureFeatureSelection( pdf_classes, fast_search ) ;  
% %         id_cost = 2 ;
% %         rslt = Cost(3:size(Cost,1),:) ;
% % %         
% % %         max_v = (rslt(id_cost,size(rslt,2))) ;
% % %         
% % %         
 
        d_b = sort(d_before)  ;
        d_n = sort(sel_feats)  ;
        if isequal(d_b, d_n) || isempty(d_before)
            
        else            
            if ~isequal(d_b, d_n) %length(d_b) < length(d_n)
                hyper_output_kde_cl.sub_feature_sel_forgetting = 20*length(hyper_input_kde_cl.kde_cl) ;
            else
                hyper_output_kde_cl.sub_feature_sel_forgetting = 0 ;
            end
        end        
    case 'auto_choose_model'
        % get confusion matrix
        rslt = executeOperatorIKDEClsfr( hyper_input_kde_cl, 'introspect' ) ;
        
        % create distribution
%         C_m = rslt.Con_matrix(1:size(rslt.Con_matrix(:,1),1),1:size(rslt.Con_matrix(:,1))) ;
%         P = (1 - diag(C_m))+1e-20 ;          
%         P = exp((1 - rslt.Conf_array).^2) ;
 
        P = rslt.Conf_array.^2 ;
        for i = 1 : length(exclude_model)
            P(str2double(exclude_model{i})) = 0 ;
        end
        
        P = P / max(sum(P),1e-30) ;
                
        
        % sample class from the distribution
%         if rand(1) > 0.1
            cls = sampleDiscrete( P ) ;
%         else
%             [valp,cls] = max(P) ;
 
        rslt = executeOperatorIKDEClsfr( hyper_input_kde_cl, 'get_name_at_index', cls ) ;
        hyper_output_kde_cl = [] ;
        hyper_output_kde_cl.class = cls ;
        hyper_output_kde_cl.class_name = rslt ;                
    case 'autonomous_update' ;
        answers = [] ;
 
        for i = 1 : length(input_data) 
            Ccost = 0 ;
            rslt = executeOperatorIKDEClsfr( hyper_input_kde_cl, 'input_data', input_data{i}.data, 'classifyData',...
                                             'use_unknown_model',1, 'extensive_answer', 1  ) ;
            if isequal( autonomous_update, 'oracle_verified' )
                % check if the label already exists and check for errors
                [data, class, class_name, class_exists] = parseClassData(hyper_input_kde_cl, input_data{i} ) ;
                answers = horzcat(answers, 0) ;
                if rslt.H >  hyper_input_kde_cl.autoUpdateThres.upper || rslt.C ~= class                    
                    hyper_input_kde_cl = executeOperatorIKDEClsfr( hyper_input_kde_cl, 'input_data', input_data(i), 'add_input', vforwvargin{:} ) ;
                end
            elseif isequal( autonomous_update, 'self_verified' )
                if rslt.H >  hyper_input_kde_cl.autoUpdateThres.upper || rslt.C == -1 % uncertain classification 
                    hyper_input_kde_cl = executeOperatorIKDEClsfr( hyper_input_kde_cl, 'input_data', input_data(i), 'add_input', vforwvargin{:} ) ;
                    answers = horzcat(answers, Cdescribe) ;
                else
                    %Certain Classification 
                    if rslt.H < hyper_input_kde_cl.autoUpdateThres.lower
                        answers = horzcat(answers, 0) ;
                        input_data{i}.class_name = [] ;
                        input_data{i}.class = rslt.C ; 
                        hyper_input_kde_cl = executeOperatorIKDEClsfr( hyper_input_kde_cl, 'input_data', input_data(i), 'add_input', vforwvargin{:} ) ;
                    else
                        answers = horzcat(answers, 0) ;
                    end
                end
            elseif isequal( autonomous_update, 'situated_verified' )
                if rslt.H > 1e-6 || rslt.C == -1 %rslt.H >  hyper_input_kde_cl.autoUpdateThres.upper || rslt.C == -1 % uncertain classification 
                    Ccost = Ccost + Clist ;
                    if rslt.C == -1  || max(rslt.P) < 0.2  % ask open question
%                         hyper_input_kde_cl = executeOperatorIKDEClsfr( hyper_input_kde_cl, 'input_data', input_data(i), 'add_input', vforwvargin{:} ) ;
                        hyper_input_kde_cl = unlearn_few_update_correct(hyper_input_kde_cl, rslt,...
                                                                         input_data(i), 0.05, vforwvargin) ;    
       
                        Ccost = Ccost + Cdescribe ;
                    else % ask polar question 
                        % is this class rslt.C?
                        Ccost = Ccost + Cans ;
                        Cl_name = executeOperatorIKDEClsfr( hyper_input_kde_cl, 'get_name_at_index', rslt.C) ;
                        if isequal(input_data{i}.class_name, Cl_name) % yes
%                             hyper_input_kde_cl = executeOperatorIKDEClsfr( hyper_input_kde_cl, 'input_data', input_data(i), 'add_input', vforwvargin{:} ) ;                            
                                hyper_input_kde_cl = unlearn_few_update_correct(hyper_input_kde_cl, rslt,...
                                                                         input_data(i), 0.05, vforwvargin) ;
                        else % no           
                              hyper_input_kde_cl = unlearn_few_update_correct(hyper_input_kde_cl, rslt,...
                                                                         input_data(i), 0.05, vforwvargin) ; 
%                             cls_tmp = input_data{i}.class_name ;
%                             input_data{i}.class_name = Cl_name ; 
%                             hyper_input_kde_cl = executeOperatorIKDEClsfr( hyper_input_kde_cl, 'input_data', input_data(i), 'unlearn_with_input', vforwvargin{:} ) ;
%                             input_data{i}.class_name = cls_tmp ;
%                             hyper_input_kde_cl = executeOperatorIKDEClsfr( hyper_input_kde_cl, 'input_data', input_data(i), 'add_input', vforwvargin{:} ) ;
                            Ccost = Ccost + Cdescribe ;
                        end                        
                    end
%                     Ccost = 1 ;
                    answers = horzcat(answers, Ccost) ;
                else % Certain Classification so no questions
%                    if rslt.H <= hyper_input_kde_cl.autoUpdateThres.lower
%                         answers = horzcat(answers, 0) ;
%                         input_data{i}.class_name = [] ;
%                         input_data{i}.class = rslt.C ; 
%                         hyper_input_kde_cl = executeOperatorIKDEClsfr( hyper_input_kde_cl, 'input_data', input_data(i), 'add_input', vforwvargin{:} ) ;
%                    else
                       answers = horzcat(answers, 0) ;
%                    end
                end
            else
                error('Unknown "autonomous_update" type!') ;
            end
        end       
        hyper_output_kde_cl = hyper_input_kde_cl ;
        hyper_output_kde_cl.answers = answers ;
    case 'add_input'
        % determine if the class already exists, and initialize it if it
        % does not exist

% % % % %          
% % % % % if length(hyper_input_kde_cl.kde_cl) > 0
% % % % %             pdf.Mu = [] ;
% % % % %             pdf.Cov = {} ;
% % % % %             pdf.w = [] ;
% % % % %             N_eff = 0 ;
% % % % %             for i = 1 : length(hyper_input_kde_cl.kde_cl) 
% % % % %                pdf.Mu = horzcat(pdf.Mu, hyper_input_kde_cl.kde_cl{i}.pdf.Mu) ;
% % % % %                pdf.Cov = horzcat(pdf.Cov, hyper_input_kde_cl.kde_cl{i}.pdf.smod.ps.Cov) ; 
% % % % %                pdf.w = horzcat(pdf.w, hyper_input_kde_cl.kde_cl{i}.pdf.w) ;  
% % % % %                N_eff = N_eff + hyper_input_kde_cl.kde_cl{i}.ikdeParams.N_eff ;
% % % % %             end
% % % % %             [new_mu, new_Cov, new_w] = momentMatchPdf(pdf.Mu, pdf.Cov, pdf.w) ;
% % % % %             
% % % % %             H = ndDirectPlugin_JointClean( pdf.Mu, pdf.Cov, pdf.w, new_Cov, N_eff ) ;
% % % % %             hyper_input_kde_cl.kde_cl{i}.pdf.smod.H = H ;
% % % % % end
        
        
        answers = [] ;
        for i = 1 : length(input_data)  
            if isempty(input_data{i}.data)
                continue ; 
            end
            Ccost = 0 ; 
            % Ccost = 1 ; 
            answers = horzcat(answers, Ccost*ones(1,size(input_data{i}.data,2))) ;
            % check if the label already exists and check for errors
            [data, class, class_name, class_exists] = parseClassData(hyper_input_kde_cl, input_data{i} ) ;
            
            % initialize cumulative feature selector
            if isempty(hyper_input_kde_cl.cummulative_feat_costs)
               hyper_input_kde_cl.cummulative_feat_costs = zeros(1,size(data,1)) ;  
            end
            
            % if attenuation should be triggered
            if hyper_input_kde_cl.sub_feature_sel_forgetting > 0
                hyper_input_kde_cl.sub_feature_sel_forgetting = hyper_input_kde_cl.sub_feature_sel_forgetting - 1 ;
                kde_w_attenuation = 1 - 1/10  ;
            else
                kde_w_attenuation = 1  ;
            end
            kde_w_attenuation = 1 ;
 
            if class_exists == 0                                        
                % initialize new class    
                kde = executeOperatorIKDE( [], 'input_data', input_data{i}.data, 'add_input', ...                                           
                                           'compressionClusterThresh', hyper_input_kde_cl.compressionClusterThresh,...
                                           'maxNumCompsBeforeCompression', hyper_input_kde_cl.force_value_init_of_maxNumCompsBeforeCompression,...
                                           'turn_off_splitting', hyper_input_kde_cl.turn_off_splitting,...
                                           vforwvargin{:}) ; 
  
                hyper_input_kde_cl.kde_cl = horzcat(hyper_input_kde_cl.kde_cl, kde) ;                
                hyper_input_kde_cl.class_labels = horzcat(hyper_input_kde_cl.class_labels, class) ;
                hyper_input_kde_cl.class_labels_names = horzcat(hyper_input_kde_cl.class_labels_names, class_name) ;
                hyper_input_kde_cl.unknown_model_value = prob_of_unknown( hyper_input_kde_cl.kde_cl ) ;
            else
                % update existing class
            
% %                 % if minimal_required_examps_mode is on, then first check
% %                 % for classification
% %                 if hyper_input_kde_cl.minimal_required_examps_mode == 1
% %                     reslt = executeOperatorIKDEClsfr( hyper_input_kde_cl, 'input_data', data, 'classifyData', 'use_unknown_model', 1 ) ;                     
% %                     if reslt.C == class
% %                        % correct classification, so continue  
% %                        hyper_input_kde_cl.kde_cl{class}.ikdeParams.N_eff = hyper_input_kde_cl.kde_cl{class}.ikdeParams.N_eff + 1;
% %                        hyper_input_kde_cl.kde_cl{class} = executeOperatorIKDE( hyper_input_kde_cl.kde_cl{class}, 'recalculate_bandwidth', ...
% %                                                                                 'input_data', [], 'add_input') ;                       
% %                        continue ; 
% %                     end                    
% %                 end               

                if isequal(hyper_input_kde_cl.typeRecDescr,'dKDE')
                    % make negative data model
                    otherClasses = makeOtherClasses( hyper_input_kde_cl.kde_cl, class, size(data,2), use_equalimportance, hyper_input_kde_cl.pair_dist_struct ) ; 
                elseif isequal(hyper_input_kde_cl.typeRecDescr,'oKDE')
                    otherClasses = {} ;
                elseif isequal(hyper_input_kde_cl.typeRecDescr,'AM')
                    otherClasses = {} ;
                elseif isequal(hyper_input_kde_cl.typeRecDescr,'dAM')
                    otherClasses = makeOtherClasses( hyper_input_kde_cl.kde_cl, class, size(data,2), use_equalimportance, hyper_input_kde_cl.pair_dist_struct ) ; %input_kde_cl  output_kde_cl
                else
                    error('Unknown update rule! Either reconstructive or discriminative !') ;
                end
                
                sub_feats = [] ;
                if hyper_input_kde_cl.react_compression_to_feature_selection == 1 
%                     &&...
%                          hyper_input_kde_cl.Params.minNumDataPointsToFormKDE + 1 > hyper_input_kde_cl.kde_cl{class}.ikdeParams.N_eff
                    sub_feats = hyper_input_kde_cl.sub_selected_features ;
                end
  
 try
     
                                                        
  
     
%                 msg = sprintf('Processing kde no: %d', class) ; disp(msg)
                % update the kde
                hyper_input_kde_cl.kde_cl{class} = ...
                                  executeOperatorIKDE( hyper_input_kde_cl.kde_cl{class}, 'input_data', ...
                                                       data, 'add_input', 'otherClasses', otherClasses,...
                                                       'selectSubDimensions', sub_feats,...
                                                       'kde_w_attenuation', kde_w_attenuation, ...
                                                       'turn_off_splitting', hyper_input_kde_cl.turn_off_splitting,...
                                                       vforwvargin{:},...
                                                       'compressionClusterThresh', hyper_input_kde_cl.compressionClusterThresh ) ;  
                                             
%                     if sum(sum(imag( hyper_input_kde_cl.kde_cl{class}.pdf.smod.H))) > 0
%             error('ERROR!!!') ;
%             
%         end                                  
                                                 
                                                   
 catch
    dfg 
 end
            end            
        end

        % search for degenerate kdes and reapproximate their bandwidths
%         for i = 1 :0 %length(hyper_input_kde_cl.kde_cl)
%            is_kde_degenerate = executeOperatorIKDE( hyper_input_kde_cl.kde_cl{i}, 'test_if_kde_is_degenerated') ; 
%            if is_kde_degenerate == 1 %hyper_input_kde_cl.kde_cl{i}.ikdeParams.N_eff < hyper_input_kde_cl.Params.minNumDataPointsToFormKDE           
%                if isequal(hyper_input_kde_cl.typeRecDescr,'dKDE') || isequal(hyper_input_kde_cl.typeRecDescr,'oKDE')
%                     otherClasses = makeOtherClasses( hyper_input_kde_cl.kde_cl, i, 0, use_equalimportance, hyper_input_kde_cl.pair_dist_struct ) ; 
%                else
%                    otherClasses = {} ;
%                end
%                
%                if isequal(hyper_input_kde_cl.typeRecDescr,'AM') 
%                    sz = size(hyper_input_kde_cl.kde_cl{i}.pdf.Cov{1}) ;
%                    for kk = 1 : length(hyper_input_kde_cl.kde_cl{i}.pdf.w)
%                        hyper_input_kde_cl.kde_cl{i}.pdf.Cov{kk} = hyper_input_kde_cl.kde_cl{i}.pdf.Cov{kk} + eye(sz)*(1e-10) ;
%                    end
%                
%                else
%                
%                if isequal(hyper_input_kde_cl.typeRecDescr,'dKDE') || isequal(hyper_input_kde_cl.typeRecDescr,'oKDE')               
%                     hyper_input_kde_cl.kde_cl{i} = ...
%                                      executeOperatorIKDE( hyper_input_kde_cl.kde_cl{i}, 'set_auxiliary_bandwidth' ,...
%                                      'otherClasses', otherClasses, ...
%                                      'turn_off_splitting', hyper_input_kde_cl.turn_off_splitting) ;
% %                else
% %                    deltaH = eye(size(hyper_input_kde_cl.kde_cl{i}.pdf.Cov{1}))*(1e-5) ;
% %                    for k = 1 : length(hyper_input_kde_cl.kde_cl{i}.pdf.w)
% %                       hyper_input_kde_cl.kde_cl{i}.pdf.Cov{k} = hyper_input_kde_cl.kde_cl{i}.pdf.Cov{k} + deltaH ;
% %                    end
%                end
% 
%                end
%            end
%             
%         end
        hyper_input_kde_cl.answers = answers ;
        hyper_output_kde_cl = hyper_input_kde_cl ;   
                    
        
%         [U,S,V] = svd(hyper_input_kde_cl.kde_cl{6}.pdf.Cov{1}) ; 
%         if sum(diag(S)') > 1e3
%             dfdfg = 56
%         end
%         
        
%         if isempty(hyper_input_kde_cl.kde_cl)
%             output_kde_cl = {} ; 
%             for i = 1 : length(input_data)
%                 %         class = input_data{i}.class ;
%                 kde = executeOperatorIKDE( [], 'input_data', input_data{i}.data, 'add_input', vforwvargin{:}  ) ;
% 
%                 output_kde_cl = horzcat(output_kde_cl, kde) ;
%                 % store class labels
%                 hyper_input_kde_cl.class_labels = horzcat(hyper_input_kde_cl.class_labels, input_data{i}.class) ; 
%             end    
%             hyper_output_kde_cl = hyper_input_kde_cl ;  
%             hyper_output_kde_cl.kde_cl = output_kde_cl ;
%             hyper_output_kde_cl.unknown_model_value = prob_of_unknown( hyper_output_kde_cl.kde_cl ) ;
%             return ;
%         else
%             % store previous
%             hyper_output_kde_cl = hyper_input_kde_cl ;
%             input_kde_cl = hyper_input_kde_cl.kde_cl ;
%             output_kde_cl = input_kde_cl ;
%                      
%                 % add input to kdes
%                 for i = 1 : length(input_data)
%                     class = input_data{i}.class ;
%                     data = input_data{i}.data ;
%                     
%                     if isempty(data)
%                         continue ;
%                     end
%                     
%                     if isequal(hyper_input_kde_cl.typeRecDescr,'dKDE')
%                         % make negative data model
%                         otherClasses = makeOtherClasses( input_kde_cl, class, size(data,2), use_equalimportance ) ; %input_kde_cl  output_kde_cl
%                     elseif isequal(hyper_input_kde_cl.typeRecDescr,'oKDE')
%                         otherClasses = {} ;
%                     elseif isequal(hyper_input_kde_cl.typeRecDescr,'AM')
%                         otherClasses = {} ;    
%                     elseif isequal(hyper_input_kde_cl.typeRecDescr,'dAM')
%                         otherClasses = makeOtherClasses( input_kde_cl, class, size(data,2), use_equalimportance ) ; %input_kde_cl  output_kde_cl
%                     else
%                         error('Unknown update rule! Either reconstructive or discriminative !') ;
%                     end
%                     
%                     % update the kde
%                     kde = executeOperatorIKDE( input_kde_cl{class}, 'input_data', ...
%                         data, 'add_input',...
%                         'otherClasses', otherClasses, vforwvargin{:} ) ;
%                     output_kde_cl{class} = kde ;
%                 end
%             hyper_output_kde_cl.kde_cl = output_kde_cl ;
%         end
    case 'unlearn_with_input'         
        for i = 1 : length(input_data)
            % check if the label already exists and check for errors
            [data, class, class_name, class_exists] = parseClassData(hyper_input_kde_cl, input_data{i} ) ;
            
%             class = input_data{i}.class ;
%             data = input_data{i}.data ;     
            
            % create negative classes
            otherClasses = {} ;
            if isequal(hyper_input_kde_cl.typeRecDescr,'dKDE')
                % make negative data model
                otherClasses = makeOtherClasses( hyper_input_kde_cl.kde_cl,...
                               class, size(data,2), use_equalimportance, hyper_input_kde_cl.pair_dist_struct ) ; %input_kde_cl              
            end

           sub_feats = [] ;
                if hyper_input_kde_cl.react_compression_to_feature_selection == 1 %&&...
%                          hyper_input_kde_cl.Params.minNumDataPointsToFormKDE + 1 > hyper_input_kde_cl.kde_cl{class}.ikdeParams.N_eff
                    sub_feats = hyper_input_kde_cl.sub_selected_features ;
                end
            
            
            % unlearn the kde
            hyper_input_kde_cl.kde_cl{class} =...
                        executeOperatorIKDE( hyper_input_kde_cl.kde_cl{class}, 'input_data', ...
                        data, 'unlearn_with_input', 'otherClasses', otherClasses, 'selectSubDimensions', sub_feats, vforwvargin{:} ) ;
%             hyper_input_kde_cl.kde_cl{class}.pdf.smod.useVbw        
        end        
        hyper_output_kde_cl = hyper_input_kde_cl ; 
    case 'compress_pdf'
        % store previous
        hyper_output_kde_cl = hyper_input_kde_cl ;
        input_kde_cl = hyper_input_kde_cl.kde_cl ;
        output_kde_cl = input_kde_cl ;
        for class = 1 : length(input_kde_cl)
            otherClasses = {} ;
            if isequal(hyper_input_kde_cl.typeRecDescr,'discriminative')
                % make negative data model
                otherClasses = makeOtherClasses( output_kde_cl, class, 0, use_equalimportance, hyper_input_kde_cl.pair_dist_struct ) ; %input_kde_cl
            end
            % update the kde
            kde = executeOperatorIKDE( input_kde_cl{class}, 'compress_pdf',...
                'otherClasses', otherClasses, vforwvargin{:} ) ;
            output_kde_cl{class} = kde ;
        end
        hyper_output_kde_cl.kde_cl = output_kde_cl ;
    case 'classifyData' % classifies in the common subspace
        input_kde_cl = hyper_input_kde_cl.kde_cl ;
                P = zeros(length(input_kde_cl),size(input_data,2)) ;
                p_pr = zeros(1,length(input_kde_cl)) ;
        %         allClasses = makeOtherClasses( input_kde_cl, i, 0, use_equalimportance, tmp_pd_Struct  ) ;
        
        
        % calculate Ht from all classes
        pdf_tmp.Mu = [] ;
        pdf_tmp.Cov = {} ;
        pdf_tmp.w = [] ;
        N_eff = 0 ;
        for i = 1 : length(input_kde_cl)
            pdf_tmp.Mu = horzcat(pdf_tmp.Mu,input_kde_cl{i}.ikdeParams.scale.Mu) ;
            pdf_tmp.Cov = horzcat(pdf_tmp.Cov,input_kde_cl{i}.ikdeParams.scale.Cov) ;
            pdf_tmp.w = horzcat(pdf_tmp.w, 1) ;
            N_eff = N_eff + input_kde_cl{i}.ikdeParams.N_eff ;
        end
        pdf_tmp.w = pdf_tmp.w / sum(pdf_tmp.w) ;
        [new_mu, new_Cov, w_out] = momentMatchPdf(pdf_tmp.Mu, pdf_tmp.Cov, pdf_tmp.w) ;
        d = size(pdf_tmp.Mu,1) ;
        Ht_global = 0.2^2*(new_Cov *(4/((d+2)*N_eff))^(2/(d+4))) ;
        
        
        [U,S,V] = svd(new_Cov) ;
        V = U ;
%         id_valid = 1 : length(diag(S)) ;
%         id_valid = diag(S) >  1e-10 ;
        s = diag(S) ;
        id_valid = s/max(s) >  1e-7; % 0; -1 ; %%
        S = S(id_valid, id_valid) ;
        V = V(:,id_valid) ;
        
        % transform into common subspace
        F_trns = sqrt(inv(S))* V' ; %eye(size(new_Cov)); %
        for i = 1 : length(input_kde_cl)
%             tmp_pd_Struct = hyper_input_kde_cl.pair_dist_struct ;
%                 tmp_pd_Struct.use_approx = 0 ;
%              otherClasses = makeOtherClasses( input_kde_cl, i, 0, use_equalimportance, tmp_pd_Struct  ) ;
%             input_kde_cl{i} = executeOperatorIKDE( input_kde_cl{i}, 'set_auxiliary_bandwidth' ,  'otherClasses', otherClasses) ;
            
            input_kde_cl{i}.pdf.Mu = F_trns*(input_kde_cl{i}.pdf.Mu - repmat(new_mu,1,length(input_kde_cl{i}.pdf.w))) ;
            
              
%                 Cb = input_kde_cl{i}.pdf.smod.H + eye(size(input_kde_cl{i}.pdf.smod.H))*max(diag(input_kde_cl{i}.pdf.smod.H))*1e-1 ;
                
                
%                   [U,S,V] = svd(input_kde_cl{i}.pdf.smod.H) ;
%                   s = diag(S) ; id = s /max(s) < 1e-2 ;
%                   s(id) = max(diag(input_kde_cl{i}.pdf.smod.H))*1e-1 ;
%                  Cb = input_kde_cl{i}.pdf.smod.H  ;  
                  
                  [U,S,V] = svd(input_kde_cl{i}.pdf.smod.H) ;
                  Sg = U'*Ht_global*V ;                  
                  s = diag(S) ; id = s /max(s) < 1e-7 ;
                  sg = diag(Sg) ; s(id) = s(id) + sg(id) ;
                  S = diag(s) ;
                  input_kde_cl{i}.pdf.smod.H = U*S*V' ;
%                   Cb = input_kde_cl{i}.pdf.smod.H  ;
 
 
               [U,S,V] = svd(input_kde_cl{i}.pdf.smod.H);
                s = diag(S) ; id_invalid = s/max(s) < 1e-7 ; %1e-10 ;
                s(id_invalid) = max(s)*(1e-7) ; 
                Cb = U*diag(s)*U' ;
                  
            
            for j = 1 : length(input_kde_cl{i}.pdf.w)
                
       
%                 [U,S,V] = svd(input_kde_cl{i}.pdf.smod.H + eye(size(input_kde_cl{i}.pdf.smod.H))*max(diag(input_kde_cl{i}.pdf.smod.H))*1e-7) ; 
%                 S = U'*input_kde_cl{i}.pdf.smod.H*V ;
                
%                 [U,S,V] = svd(input_kde_cl{i}.pdf.smod.H);
%                 s = diag(S) ; id_invalid = s/max(s) < 1e-7 ; %1e-10 ;
%                 s(id_invalid) = max(s)*(1e-7) ; 
%                 Cb = U*diag(s)*U' ;
                




%                 [Utt, Stt, Vtt] = svd(input_kde_cl{i}.pdf.Cov{j}) ; Pt = diag(Stt)' ;
            
                input_kde_cl{i}.pdf.Cov{j} = input_kde_cl{i}.pdf.smod.ps.Cov{j} + Cb ;
                
               
%              [Ut, St, Vt] = svd(input_kde_cl{i}.pdf.Cov{j}) ; Pt = [Pt; diag(St)' ]
   
                
                
                input_kde_cl{i}.pdf.Cov{j} = F_trns*input_kde_cl{i}.pdf.Cov{j}*F_trns' ;
   
                
            end
            input_kde_cl{i}.pdf.Mu = input_kde_cl{i}.pdf.Mu(id_valid,:) ;
        end
        input_data = F_trns*(input_data  - repmat(new_mu,1,size(input_data,2)));
        input_data = input_data(id_valid,:) ;
        
        
        for i = 1 : length(input_kde_cl)
            p_pr(i) = input_kde_cl{i}.ikdeParams.N_eff ;
            reslt = evaluatePointsUnderPdf( input_kde_cl{i}.pdf, input_data ) ;
            P(i,:) = reslt ;
        end
        class_labels = hyper_input_kde_cl.class_labels ;
        p_pr = p_pr / sum(p_pr) ;
        if classifyWithPosterior == 1
            for i = 1 : length(input_kde_cl)
                P(i,:) = P(i,:) * p_pr(i) ;
            end
        end
        
        
        hyper_output_kde_cl = [] ;
        hyper_output_kde_cl.P = P ;
        % if we are using the unknown model then the set of classes is
        % augmented by an "unknown label" -1 !
        class_labels = hyper_input_kde_cl.class_labels ;
        if hyper_input_kde_cl.use_unknown_model == 1
            hyper_input_kde_cl.unknown_model_value = prob_of_unknown( hyper_input_kde_cl.kde_cl ) ;
            
            Rc = hyper_input_kde_cl.unknown_model_value*ones(1,size(P,2)) ;%   repmat(hyper_input_kde_cl.unknown_model_value, size(P,1), 1) ;
            for i = 1 : size(P,1)
               P(i,:) = P(i,:).*(1-Rc) ; 
            end
            
            P = [P; Rc] ;
            class_labels = horzcat(hyper_input_kde_cl.class_labels, -1) ;
        end
        
        
        sP = sum(P) ;
        for i = 1 : size(P,1)
            P(i,:) = P(i,:) ./ sP ;
        end
        
        hyper_output_kde_cl.P = P ;
        hyper_output_kde_cl.sP = sP ;
        [vals, i_max] = max(P) ;
        hyper_output_kde_cl.C = class_labels(i_max) ;
        
        if extensive_answer == 1
           % calculate and output additional statistics
           num_classes = length(class_labels) ;
           % 1. get entropy
           H = 2*(1-sum(P.*P))  / (2*(1-1/num_classes));
           
           % 2. get R1 and R2 score (likelihood ratio)
           srtP = sort(P,1, 'descend') ;
           R1 = srtP(1,:)./ (sum(srtP(2:size(srtP,1),:)) + 1e-30) ;
           R2 = srtP(1,:)./ (srtP(2,:)+1e-30)  ;
           
           hyper_output_kde_cl.H = H ;
           hyper_output_kde_cl.R1 = R1 ;
           hyper_output_kde_cl.R2 = R2 ;
        end
        
    case 'classifyDataPreviousVersion'
        input_kde_cl = hyper_input_kde_cl.kde_cl ;       
        P = zeros(length(input_kde_cl),size(input_data,2)) ;
        p_pr = zeros(1,length(input_kde_cl)) ;
        for i = 1 : length(input_kde_cl)
            if sum(i == ignoreClasses) > 0
               continue ; 
            end
            p_pr(i) = input_kde_cl{i}.ikdeParams.N_eff ;
            
            sub_feats = [] ;
            if hyper_input_kde_cl.react_compression_to_feature_selection == 1
                sub_feats = hyper_input_kde_cl.sub_selected_features ;
            end 
            sub_feats = hyper_input_kde_cl.sub_selected_features ;
            
            if isequal(hyper_input_kde_cl.typeRecDescr,'AM') || isequal(hyper_input_kde_cl.typeRecDescr,'dAM')
                pdf_t = input_kde_cl{i}.pdf ;
%                 dd = size(pdf_t.Cov{1},1) ;
%                 for ik = 1 : length(pdf_t.Cov)
%                     pdf_t.Cov{ik} = pdf_t.Cov{ik}+eye(dd)*(1e-3)^2 ;                           
%                 end


%                 [pdf_in, data, svdRes] = regulrMixtureIntoSubspace(pdf_t, input_data) ;
%                 [pdf_in, data, svdRes] = projectMixtureIntoSubspace(pdf_t, input_data ) ;
%                 if ~isempty(svdRes)
%                     reslt.evalpdf = evaluatePointsUnderPdf(pdf_in, data) ;        
%                 else
                    reslt.evalpdf = evaluatePointsUnderPdf(pdf_t, input_data) ;
%                 end
            else    
                
                % search for degenerate kdes and reapproximate their bandwidths
                is_kde_degenerate = executeOperatorIKDE( hyper_input_kde_cl.kde_cl{i}, 'test_if_kde_is_degenerated') ;
                tmp_pd_Struct = hyper_input_kde_cl.pair_dist_struct ;
                tmp_pd_Struct.use_approx = 0 ;
                if is_kde_degenerate == 1 %hyper_input_kde_cl.kde_cl{val_get}.ikdeParams.N_eff < hyper_input_kde_cl.Params.minNumDataPointsToFormKDE
                    otherClasses = makeOtherClasses( input_kde_cl, i, 0, use_equalimportance, tmp_pd_Struct  ) ;
                    input_kde_cl{i} = executeOperatorIKDE( input_kde_cl{i}, 'set_auxiliary_bandwidth' ,  'otherClasses', otherClasses) ;
                end
                reslt = executeOperatorIKDE( input_kde_cl{i}, 'input_data',  input_data, 'evalPdfOnData',...
                    'selectSubDimensions', sub_feats) ;                
            end
            P(i,:) = reslt.evalpdf ; %p ;
        end
        p_pr = p_pr / sum(p_pr) ;
        if classifyWithPosterior == 1
            for i = 1 : length(input_kde_cl)
                P(i,:) = P(i,:) * p_pr(i) ;
            end            
        end
    
        
        hyper_output_kde_cl = [] ;
        hyper_output_kde_cl.P = P ;
  
        % if we are using the unknown model then the set of classes is
        % augmented by an "unknown label" -1 !
        class_labels = hyper_input_kde_cl.class_labels ;
        if hyper_input_kde_cl.use_unknown_model == 1
            hyper_input_kde_cl.unknown_model_value = prob_of_unknown( hyper_input_kde_cl.kde_cl ) ;
            
            Rc = hyper_input_kde_cl.unknown_model_value*ones(1,size(P,2)) ;%   repmat(hyper_input_kde_cl.unknown_model_value, size(P,1), 1) ;
            for i = 1 : size(P,1)
               P(i,:) = P(i,:).*(1-Rc) ; 
            end
            
            P = [P; Rc] ;
            class_labels = horzcat(hyper_input_kde_cl.class_labels, -1) ;
        end
        
        sP = sum(P) ;
        for i = 1 : size(P,1)
            P(i,:) = P(i,:) ./ sP ;
        end
        
        hyper_output_kde_cl.P = P ;
        hyper_output_kde_cl.sP = sP ;
        [vals, i_max] = max(P) ;        
        hyper_output_kde_cl.C = class_labels(i_max) ;
        
        if extensive_answer == 1
           % calculate and output additional statistics
           num_classes = length(class_labels) ;
           % 1. get entropy
           H = 2*(1-sum(P.*P))  / (2*(1-1/num_classes));
           
           % 2. get R1 and R2 score (likelihood ratio)
           srtP = sort(P,1, 'descend') ;
           R1 = srtP(1,:)./ (sum(srtP(2:size(srtP,1),:)) + 1e-30) ;
           R2 = srtP(1,:)./ (srtP(2,:)+1e-30)  ;
           
           hyper_output_kde_cl.H = H ;
           hyper_output_kde_cl.R1 = R1 ;
           hyper_output_kde_cl.R2 = R2 ;
        end
    case 'calculate_gains'
        % calculate gains for all input data-points
        hyper_output_kde_cl = {} ;
        for i = 1 : size(input_data,2)        
            % classify data_point
            rslt = executeOperatorIKDEClsfr( hyper_input_kde_cl, 'input_data', input_data(:,i), ...
                                             'classifyData', 'use_unknown_model', 1,...
                                             'extensive_answer', 1 ) ;
            
            I1 = rslt.H ; % .R1, .R2
            P1 = rslt.P ;
            C1 = rslt.C(1) ; 
            
            I2 = 0 ; I3 = 0 ;
            % after first update
            for j = 1 : length(hyper_input_kde_cl.kde_cl)              
                dat.class = j ;
                dat.data = input_data(:,i) ;           
                kde_cl_tmp = executeOperatorIKDEClsfr( hyper_input_kde_cl, 'input_data', {dat}, 'add_input', 'force_prevent_compression', 1 ) ;
                
                rslt = executeOperatorIKDEClsfr( kde_cl_tmp, 'input_data', input_data(:,i), ...
                                             'classifyData', 'use_unknown_model', 1,...
                                             'extensive_answer', 1 ) ;
                dI2 = rslt.H ;
                I2 = I2 + dI2*P1(j) ;
                P2 = rslt.P ;
                
                dI3 = 0 ;
                for k = 1 : length(hyper_input_kde_cl.kde_cl)
                    dat.class = k ;
                    dat.data = input_data(:,i) ;           
                    kde_cl_tmp2 = executeOperatorIKDEClsfr( kde_cl_tmp, 'input_data', {dat}, 'add_input', 'force_prevent_compression', 1 ) ;
                    rslt = executeOperatorIKDEClsfr( kde_cl_tmp2, 'input_data', input_data(:,i), ...
                                                    'classifyData', 'use_unknown_model', 1,...
                                                    'extensive_answer', 1 ) ;
                    ddI3 = rslt.H ;
                    dI3 = dI3 + ddI3*P2(k) ;
                end
                dI3 = dI3 / sum(P2(1:length(P2)-1)) ;
                I3 = I3 + dI3*P1(j) ;
            end
            I2 = I2 / sum(P1(1:length(P1)-1)) ; 
            I3 = I3 / sum(P1(1:length(P1)-1)) ; 
            
            rslt = [] ;
            rslt.I1 = I1 ;
            rslt.I2 = I2 ;
            rslt.I3 = I3 ;
            rslt.Gain1 = I1 - I2 ;
            rslt.Gain2 = I2 - I3 ;
            rslt.dGain = rslt.Gain2 - rslt.Gain1 ;
            rslt.C = C1 ;
            hyper_output_kde_cl = horzcat(hyper_output_kde_cl, rslt) ;
        end 
        
    case 'predict_missing'
        input_kde_cl = hyper_input_kde_cl.kde_cl ;
        hyper_output_kde_cl = {} ;
        for i = 1 : length(input_kde_cl)
            reslt = executeOperatorIKDE( input_kde_cl{i}, 'input_data', input_data, 'predictMissingVals', vforwvargin{:}  ) ;
            hyper_output_kde_cl = horzcat(hyper_output_kde_cl, reslt) ;
        end   
    case 'introspect'
        % generate samples from each model, classify them and compute the
        % confusion matrix Con_matrix
        use_missclassification_in_gains = 1 ;
        % get auxilliary bandwidth
        w = [] ; Mu = [] ; Cov = {} ; Neff = 0 ;
        for i = 1 : length(hyper_input_kde_cl.kde_cl)
            w = [w, hyper_input_kde_cl.kde_cl{i}.ikdeParams.N_eff*hyper_input_kde_cl.kde_cl{i}.pdf.w] ;
            Mu = [Mu, hyper_input_kde_cl.kde_cl{i}.pdf.Mu] ;
            Cov = horzcat(Cov, hyper_input_kde_cl.kde_cl{i}.pdf.Cov) ;    
            Neff = Neff + hyper_input_kde_cl.kde_cl{i}.ikdeParams.N_eff ;
        end
        w = w / sum(w) ;
        
        % replace bandwidths if necessary
        d = size(Mu,1) ;
        [new_mu, new_Cov, w_out] = momentMatchPdf(Mu, Cov, w) ;
        nH = new_Cov *(4/((d+2)*Neff))^(2/(d+4)) ;
        for i = 1 : length(hyper_input_kde_cl.kde_cl)
%             if (hyper_input_kde_cl.kde_cl{i}.ikdeParams.N_eff/4 < min([d,length(hyper_input_kde_cl.sub_selected_features)])) 
            lmb = min([d,length(hyper_input_kde_cl.sub_selected_features)]) ;
                w1 = exp(- hyper_input_kde_cl.kde_cl{i}.ikdeParams.N_eff/ lmb ) ;
                w2 = 1 - w1 ;
               hyper_input_kde_cl.kde_cl{i}.pdf.smod.H = nH*w1 + hyper_input_kde_cl.kde_cl{i}.pdf.smod.H*w2 ; 
               hyper_input_kde_cl.kde_cl{i}.pdf = getKDEfromSampleDistribution( hyper_input_kde_cl.kde_cl{i}.pdf ) ;                  
%             end           
        end
  

        N_samps = 100*max([1,length(hyper_input_kde_cl.sub_selected_features)]) ; 
        idx_unknown = length(hyper_input_kde_cl.kde_cl) + 1 ;
        Con_matrix = zeros(length(hyper_input_kde_cl.kde_cl), idx_unknown) ;
        Conf_array = zeros(1, length(hyper_input_kde_cl.kde_cl)) ;
        for i = 1 : length(hyper_input_kde_cl.kde_cl)         
            x_tmp = sampleGaussianMixture( hyper_input_kde_cl.kde_cl{i}.pdf, N_samps ) ;
            rslt = executeOperatorIKDEClsfr( hyper_input_kde_cl, 'input_data', x_tmp, ...
                      'classifyData', 'use_unknown_model', 1, 'extensive_answer', 1  ) ;%'use_unknown_model', 0
            for j = 1 : length(hyper_input_kde_cl.kde_cl)                               
               Con_matrix(i,j) = sum(rslt.C == j)/N_samps ;                
            end  
            H = sort(rslt.H, 'ascend') ;            
            Conf_array(i) = mean( H ) ;%H(round(N_samps/2):length(H)) ) ; 
            Con_matrix(i,idx_unknown) = sum(rslt.C == -1)/N_samps ;
        end
        
        if use_missclassification_in_gains == 1 
            CC = Con_matrix(1:length(hyper_input_kde_cl.kde_cl), 1:length(hyper_input_kde_cl.kde_cl)) ;
            Conf_array = 1 - diag( CC )' ;
        end
        
        hyper_output_kde_cl.Conf_array = Conf_array ;
        hyper_output_kde_cl.N_samps = N_samps ;
        hyper_output_kde_cl.Con_matrix = Con_matrix ;
        
    case 'showKDE_of_class_index'        
        if val_get > length(hyper_input_kde_cl.class_labels)
            return ;
        end        
      
        sub_feats = [] ;
        if hyper_input_kde_cl.react_compression_to_feature_selection == 1
            sub_feats = hyper_input_kde_cl.sub_selected_features ;
        end
 
        if length(sub_feats) > 3
            sub_feats = sub_feats(1:3) ;                       
        end
        
        if isempty(sub_feats) && size(hyper_input_kde_cl.kde_cl{val_get}.pdf.Mu,1) > 3
            sub_feats = 1:3 ;  
        end
 
        % search for degenerate kdes and reapproximate their bandwidths
        is_kde_degenerate = executeOperatorIKDE( hyper_input_kde_cl.kde_cl{val_get}, 'test_if_kde_is_degenerated') ;
        tmp_pd_Struct = hyper_input_kde_cl.pair_dist_struct ;
        tmp_pd_Struct.use_approx = 0 ;
        if is_kde_degenerate == 1 %hyper_input_kde_cl.kde_cl{val_get}.ikdeParams.N_eff < hyper_input_kde_cl.Params.minNumDataPointsToFormKDE                  
                    otherClasses = makeOtherClasses( hyper_input_kde_cl.kde_cl, val_get, 0, use_equalimportance, tmp_pd_Struct  ) ; 
                    hyper_input_kde_cl.kde_cl{val_get} = ...
                                     executeOperatorIKDE( hyper_input_kde_cl.kde_cl{val_get}, 'set_auxiliary_bandwidth' ,...
                                     'otherClasses', otherClasses) ;
        end
        
        executeOperatorIKDE( hyper_input_kde_cl.kde_cl{val_get}, 'showKDE', 'selectSubDimensions', sub_feats, vforwvargin{:} ) ;
        hyper_output_kde_cl = hyper_input_kde_cl.class_labels_names(val_get) ;
%         title(hyper_output_kde_cl) ;  
    case 'regularize_all_kdes'
        
        % search for degenerate kdes and reapproximate their bandwidths
        tmp_pd_Struct = hyper_input_kde_cl.pair_dist_struct ;
        tmp_pd_Struct.use_approx = 0 ;
        for val_get = 1 : length(hyper_input_kde_cl.kde_cl)
            is_kde_degenerate = executeOperatorIKDE( hyper_input_kde_cl.kde_cl{val_get}, 'test_if_kde_is_degenerated') ;            
            if is_kde_degenerate == 1 %hyper_input_kde_cl.kde_cl{val_get}.ikdeParams.N_eff < hyper_input_kde_cl.Params.minNumDataPointsToFormKDE     
         
                    otherClasses = makeOtherClasses( hyper_input_kde_cl.kde_cl, val_get, 0, use_equalimportance, tmp_pd_Struct ) ; 
                    hyper_input_kde_cl.kde_cl{val_get} = ...
                                     executeOperatorIKDE( hyper_input_kde_cl.kde_cl{val_get}, 'set_auxiliary_bandwidth' ,...
                                     'otherClasses', otherClasses) ;
 
                    
            end
        end
        hyper_output_kde_cl = hyper_input_kde_cl ;
    otherwise
        hyper_output_kde_cl = hyper_input_kde_cl ;
end

% -------------------------------------------------------------------- %
function hyper_input_kde_cl = unlearn_few_update_correct(hyper_input_kde_cl, rslt, input_data, scaleth, vforwvargin)
 
Cl_index = executeOperatorIKDEClsfr( hyper_input_kde_cl, 'get_index_at_name', input_data{1}.class_name) ;
cls_tmp = input_data{1}.class_name ;
 
f = rslt.P(Cl_index) ;
try
dp = (rslt.P - f*scaleth) > 0 ; 
dp(end) = 0  ;
dp(Cl_index) = 0 ;
 
dp = find(dp) ;
catch
    dp = [] ;
end

for i = 1 : length(dp)
    input_data{1}.class_name = executeOperatorIKDEClsfr( hyper_input_kde_cl, 'get_name_at_index', dp(i)) ; 
    hyper_input_kde_cl = executeOperatorIKDEClsfr( hyper_input_kde_cl, 'input_data', input_data, 'unlearn_with_input', vforwvargin{:} ) ;    
end

input_data{1}.class_name = cls_tmp ;
hyper_input_kde_cl = executeOperatorIKDEClsfr( hyper_input_kde_cl, 'input_data', input_data, 'add_input', vforwvargin{:} ) ;
 
% ------------------------------------------------------------------ %
function otherClasses = makeOtherClasses( input_kde_cl, i_exclude, newadds, use_equalimportance, pair_dist_struct )
 
otherClasses.pdfs = {} ; 
otherClasses.N_eff = [] ;
w_other = ones(1,length(input_kde_cl)-1) ;  
w_other = w_other / sum(w_other) ;
% otherClasses.priors = []; %(length(input_kde_cl)-1) / length(input_kde_cl) ;  
otherClasses.inner_priors = [] ; %w_other ;
N_eff_all = 0 ;
N_eff_neg = 0 ;
for i = 1 : length(input_kde_cl) 
    N_eff_all = N_eff_all + input_kde_cl{i}.ikdeParams.N_eff ;
    if i ~= i_exclude
        N_eff_neg = N_eff_neg + input_kde_cl{i}.ikdeParams.N_eff ;
        if pair_dist_struct.use_approx ~= 0 
            t = test_for_overlap(input_kde_cl{i_exclude}.ikdeParams.scale.Mu, input_kde_cl{i_exclude}.ikdeParams.scale.Cov+input_kde_cl{i_exclude}.pdf.smod.H, ...
                                 input_kde_cl{i}.ikdeParams.scale.Mu, input_kde_cl{i}.ikdeParams.scale.Cov+input_kde_cl{i}.pdf.smod.H,...
                                 pair_dist_struct.dist_th) ;                             
            if t == 0
               continue ; 
            end
        end
        
        otherClasses.inner_priors = horzcat(otherClasses.inner_priors, input_kde_cl{i}.ikdeParams.N_eff ) ;
        
        pdfx = input_kde_cl{i}.pdf ; 
        pdfx.scale = input_kde_cl{i}.ikdeParams.scale ;
        otherClasses.pdfs = horzcat(otherClasses.pdfs, pdfx) ;
        otherClasses.N_eff = horzcat(otherClasses.N_eff, input_kde_cl{i}.ikdeParams.N_eff ) ;
    end
end
 
otherClasses.inner_priors = otherClasses.inner_priors / N_eff_neg ; %sum([otherClasses.inner_priors,input_kde_cl{i_exclude}.ikdeParams.N_eff + newadds ]) ;
 
otherClasses.N_all_classes = length(input_kde_cl) ;
% otherClasses.priors = sum(otherClasses.inner_priors) ;
% otherClasses.priors = 1 / length(input_kde_cl) ;
otherClasses.priors = (input_kde_cl{i_exclude}.ikdeParams.N_eff+newadds) / (N_eff_all + newadds) ;
if use_equalimportance == 1 
    otherClasses.inner_priors = w_other ;
    otherClasses.priors = (length(input_kde_cl)-1) / length(input_kde_cl) ;     
end
% if isempty(otherClasses.pdfs)
%     otherClasses = {}  ;
% end

% ------------------------------------------------------------------ %
function t = test_for_overlap(Mu_ref, Cov_ref, Mu_t, Cov_t, dist_th) 

C = Cov_ref + Cov_t + eye(size(Cov_ref))*1e-3 ;
d = Mu_ref - Mu_t ;
l = sqrt(d'*inv(C)*d) ;
t = l < dist_th ;
                             

% ------------------------------------------------------------------ %
function P = prob_of_unknown( classes_kde_cl )

N = [] ;
for i = 1 : length(classes_kde_cl)
   N = [N , classes_kde_cl{1}.ikdeParams.N_eff ] ;     
end
n = sum(N) ;
p = N/n ;

csr = n/length(classes_kde_cl) ;%P = sum(((1-p).^n).*p) ;
P = sum(((1-p).^n).*p)+0.5* +0.5*exp(-csr/10) ;
P = 1e-4 ; %max([P,1e-3]) ;

% ----------------------------------------------------------------- %
function [data, class, class_name, class_exists] = parseClassData(hyper_input_kde_cl, input_data )

class_exists = 0 ;
class = [] ;
class_name = [] ;
if isfield(input_data,'class')
    class = input_data.class ;
end
if isfield(input_data,'class_name')
    class_name = input_data.class_name ;
end
data = input_data.data ;

% continue if the data is empty
if isempty(data)
    return ;
end

% check if the label already exists and check for errors
if ~isempty(class) && class <= length(hyper_input_kde_cl.class_labels )
    class_exists = 1 ;
    if ~isempty(class_name) && ~isequal(class_name, hyper_input_kde_cl.class_labels_names{class})
        error('Error in naming the class! The consecutive number does not match the naming!') ;
    end
    class_name = hyper_input_kde_cl.class_labels_names{class} ;
end

if isempty(class)
    if ~isempty(class_name)
        for i_cls = 1 : length(hyper_input_kde_cl.class_labels_names)
            if isequal(class_name, hyper_input_kde_cl.class_labels_names{i_cls})
                class_exists = 1 ;
                class = i_cls ;
                break ;
            end
        end
    else
        error('At least one of the variables "class" or "class_name" should be specified!') ;
    end
end
% classes should be identified now

if isempty(class_name)   
    class_name = num2str(class)  ;
end
 
if isempty(class)   
    class = length(hyper_input_kde_cl.class_labels ) + 1  ;
end
 

% ---------------------------------------------------------------------- %
function [pdf_in, data, svdRes] = projectMixtureIntoSubspace(pdf_in, data )
minVals = 0.001 ;

% I = eye(size(pdf_in.Cov{1})) ;
% for i = 1 : length(pdf_in.w)
%     pdf_in.Cov{i} = pdf_in.Cov{i} + I*minVals ;    
% end
% svdRes = [] ;
% return ;
 
[new_mu, new_Cov, w_out] = momentMatchPdf(pdf_in.Mu, pdf_in.Cov, pdf_in.w) ;
[U,S,V] = svd(new_Cov) ;
V = U ;
s = diag(S) ;
 ss = s / sum(s) ;

id_valid = find(ss > minVals) ;
id_null = find(ss <= minVals) ;

if isempty(id_null)
    svdRes = [] ;
    return ;
end

id_nullVals = s(id_null) ;
d = size(new_Cov,1) ;
S_inv = eye(d,d)*0 ;
S_inv(id_valid,id_valid) = diag(diag(S(id_valid,id_valid).^(-1))) ;
for i = 1 : length(id_nullVals)
    try
        tmp_val = 1 / id_nullVals(i) ;
    catch
        tmp_val = 1 ; %1/practicallyZero ; % doesn't matter
    end
    S_inv(id_null(i),id_null(i)) = tmp_val ;
end
F_trns = sqrt(abs(S_inv))* inv(V) ;

pdf_in.Mu = F_trns*( pdf_in.Mu - repmat(new_mu,1,length(pdf_in.w))) ;
pdf_in.Mu = pdf_in.Mu(id_valid,:) ;


for i = 1 : length(pdf_in.w)
    pdf_in.Cov{i} = F_trns*pdf_in.Cov{i}*F_trns' ;
    pdf_in.Cov{i} = pdf_in.Cov{i}(id_valid,id_valid) ;
end

% forward transform the additional_data if it exists
if ~isempty(data)
    data = F_trns*(data - repmat(new_mu,1,size(data,2))) ;
    data = data(id_valid,:) ;
end

svdRes.V = V ;
svdRes.S = S ;
svdRes.id_valid = id_valid ;
svdRes.new_mu = new_mu ;
svdRes.nullspace.id_null = id_null ;
svdRes.nullspace.id_nullVals = id_nullVals ;
svdRes.nullspace.id_valid = id_valid ;
 

% ---------------------------------------------------------------------- %
function [pdf_in, data, svdRes] = regulrMixtureIntoSubspace(pdf_in, data )
minVals = 0.001 ;

% I = eye(size(pdf_in.Cov{1})) ;
% for i = 1 : length(pdf_in.w)
%     pdf_in.Cov{i} = pdf_in.Cov{i} + I*minVals ;    
% end
% svdRes = [] ;
% return ;
 
[new_mu, new_Cov, w_out] = momentMatchPdf(pdf_in.Mu, pdf_in.Cov, pdf_in.w) ;
[U,S,V] = svd(new_Cov) ;
V = U ;
s = diag(S) ;
ss = s / sum(s) ;

id_valid = find(ss > minVals) ;
id_null = find(ss <= minVals) ;

if isempty(id_null)
    svdRes = [] ;
    return ;
end

z = zeros(1, size(S,1)) ; 
z(id_null) = (1e-3)^2 ;
Z = U*diag(z)*V' ;

 
for i = 1 : length(pdf_in.w)
    pdf_in.Cov{i} =  pdf_in.Cov{i} + Z ;
end
 
svdRes.V = [] ; 
