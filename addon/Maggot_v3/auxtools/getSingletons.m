function singletons = getSingletons( pdf_ref )

singletons = [] ;
idx = [] ;
for i = 1 : length(pdf_ref.smod.q)
    if length(pdf_ref.smod.q(i).w) == 1
        idx = [idx , i] ;
    end
end

if ~isempty(idx)
    singletons = pdf_ref.Mu(:,idx) ;
end