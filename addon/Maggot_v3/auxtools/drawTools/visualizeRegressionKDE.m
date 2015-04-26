function fun_out = visualizeRegressionKDE( kde, bound, dim_a, make_spherical )

if nargin < 4
    make_spherical = 0 ;
end

n_sigma = 2 ;
n_pts = 100 ;
tabbed_x = bound(1): diff(bound)/n_pts : bound(2) ;
 
means = zeros(1, length(tabbed_x)) ;
sig_bounds = zeros(2, length(tabbed_x)) ;

for i = 1 : length(tabbed_x)
    pdf_out_approx = getConditionalKDEat( kde, tabbed_x(i), dim_a, make_spherical ) ;
    means(:,i) = pdf_out_approx.Mu ;
    
    sig_bounds(1,i) = pdf_out_approx.Mu - n_sigma*sqrt(pdf_out_approx.Cov{1}) ;
    sig_bounds(2,i) = pdf_out_approx.Mu + n_sigma*sqrt(pdf_out_approx.Cov{1}) ;   
end

x = [tabbed_x,fliplr(tabbed_x)] ;
y = [sig_bounds(1,:),fliplr(sig_bounds(2,:))] ;
 
if nargout < 1
    fun_out = [] ;
else
    fun_out = [tabbed_x ; means ] ;
end

fill(x,y,'c','LineWidth',1.5) ; hold on ;
plot(tabbed_x, means,'b','LineWidth',2) ;
