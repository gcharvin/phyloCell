function [clusters, modes_locations, modes_probs, pdf_out_moment_matched, pdf_out_l2 ] = ...
                                  detect_modes_and_approximate_by_l2_and_moment_match(pdf_ref)
% the code detects the modes on the input distribution pdf_ref
% the output is:
% clusters ... cluster indicators for the components of the input distribution 
% modes_locations ... locations of the detected modes
% modes_probs ... probability of each mode evaluated under the reference_pdf
% pdf_out_moment_matched ... compressed distribution where moment matching
%                            was used to represent a set of clustered components as a single Gaussian
% pdf_out_l2 ... compressed distribution where l2 optimization 
%                            was used to represent a set of clustered components as a single Gaussian
% Author: Matej Kristan (2011)

% detect modes by mean shift
[pdf_out_moment_matched, clusters, modes_locations] = reduceMixtureByMeanShift( pdf_ref, [], 'local_covariance' ) ;

% evaluate probability of the modes
modes_probs = evaluatePointsUnderPdf( pdf_ref, modes_locations ) ;

% construct the l2 fitted distribution
pdf_out_l2 = pdf_out_moment_matched ;
u_clusters = unique(clusters) ;
for i = 1 : length(u_clusters)
    id_selection = find(clusters==i) ;
 
    f0 = fitSingleGaussL2Approx( pdf_ref, 'idx_selection', id_selection,...
                                 'allowWeightOptimization', 1 ) ;
    pdf_out_l2.Mu(:,i) =  f0.Mu ;
    pdf_out_l2.Cov(i) =  f0.Cov ;
    pdf_out_l2.w(i) =  f0.w ;
end
pdf_out_l2.w = pdf_out_l2.w/sum(pdf_out_l2.w) ;
