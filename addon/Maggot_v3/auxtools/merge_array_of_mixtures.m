function out = merge_array_of_mixtures(input, add_pdf)

out.Mu = [] ;
out.w = [] ;
out.Cov = {} ;
for i_oth = 1 : length(input)
    out.Mu = horzcat(out.Mu, input{i_oth}.Mu) ;
    out.Cov = horzcat(out.Cov, input{i_oth}.Cov) ;
    out.w = horzcat(out.w, input{i_oth}.w) ;
end

if ~isempty(add_pdf)
    out.Mu = horzcat(out.Mu, add_pdf.Mu) ;
    out.Cov = horzcat(out.Cov, add_pdf.Cov) ;
    out.w = horzcat(out.w, add_pdf.w) ;
end
out.w = out.w/sum(out.w) ;