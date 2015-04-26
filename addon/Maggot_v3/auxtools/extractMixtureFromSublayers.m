function pdf = extractMixtureFromSublayers( kde )

pdf.Mu = [] ;
pdf.Cov = {} ;
pdf.w = [] ;
num_sublays = length(kde.smod.q) ;
for i = 1 : num_sublays
  pdf = mergeDistributions( pdf, kde.smod.q(i), [1 kde.w(i)], 0 ) ;      
end