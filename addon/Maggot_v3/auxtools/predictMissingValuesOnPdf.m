function output = predictMissingValuesOnPdf( pdf_ref, input_data, type_of_prediction )
% Predicts the missing data
% Returns:
% out.complete ... completed data  
% out.predicted ... predicted unobserved part
% out.Covincomplete ... covariance of prediction
% out.pdf ... distribution over the unobserved data

num_data = size(input_data,2) ;
output = [] ;
for i_dat = 1 : num_data
    data_i = input_data(:,i_dat) ;
    
    % 0. find missing dimensions
    id_miss = isnan(data_i) ;
    a = find(id_miss) ;
    b = find(1-id_miss) ;
    Ctol = eye(length(a))*1e-30 ;
    
    % 1. get conditional on observed values
    d = size(pdf_ref.Mu,1) ;
    len = length(pdf_ref.w) ;
    pdf_m.Mu = zeros(length(a),len) ;
    pdf_m.w = pdf_ref.w ;
    pdf_m.Cov = {} ;
    for i = 1 : len
        pdf_m.Mu(:,i) = pdf_ref.Mu(a,i) + pdf_ref.Cov{i}(a,b)*inv(pdf_ref.Cov{i}(b,b))*(data_i(b)-pdf_ref.Mu(b,i)) ;
        C = pdf_ref.Cov{i}(a,a) - pdf_ref.Cov{i}(a,b)*inv(pdf_ref.Cov{i}(b,b))*pdf_ref.Cov{i}(b,a) + Ctol ;
        pdf_m.Cov = horzcat(pdf_m.Cov, C)  ;
        
        pdf_tmp.w = 1 ;
        pdf_tmp.Mu = pdf_ref.Mu(b,i) ;
        pdf_tmp.Cov = {pdf_ref.Cov{i}(b,b)} ;
        p = evaluatePointsUnderPdf(pdf_tmp, data_i(b)) ;
        pdf_m.w(i) = pdf_m.w(i)*p ;                
    end
    pdf_m.w = pdf_m.w / sum(pdf_m.w) ;
    
    % depending on the type of predicition choose:
    switch type_of_prediction
        case 'map'
            [max_pos, max_val] = findGlobalMaximum( pdf_m ) ;
  
            out.complete = input_data(:,i_dat) ;
            out.complete(a) = max_pos ;
            out.predicted = max_pos ;
            out.Covincomplete = [] ;
            out.pdf = [] ;
        case 'expected'
            [new_mu, new_Cov, w_out] = momentMatchPdf(pdf_m.Mu, pdf_m.Cov, pdf_m.w) ;
            out.complete = input_data(:,i_dat) ;
            out.complete(a) = new_mu ;
            out.predicted = new_mu ;
            out.Covincomplete = new_Cov ;
            out.pdf = [] ;
        case 'pdf'
             out.complete = [] ; 
            out.predicted = [] ;
            out.Covincomplete = [] ;
            out.pdf = pdf_m ;            
    end
    output = horzcat(output, out) ;
end

