function I = sampleDiscrete( P )

P = P + 1e-50 ; P = P/ sum(P) ;
c = cumsum(P) ;

I = -1 ;
x = rand(1) ;
for i = 1 : length(P)
    if x <= c(i)
       I = i ;
       break ;
    end
end