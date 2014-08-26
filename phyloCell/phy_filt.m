function [ry py f]=phy_filt(x,y,mode,centerfrequency,filterwidth,filtershape,display)
% performs fft and apply custom filter to the data

% centerfrequency : frequency 
% filterwidth : width of window 
% filtershape : increasing integer make filter window more rectangular;
% default 2

% todo : apodization !

fy=fft(y);
lft1=[1:(length(fy)/2)];
lft2=[(length(fy)/2+1):length(fy)];


centerfrequency=centerfrequency*range(x);
filterwidth=filterwidth*range(x);

% Compute filter shape.
if strcmp(mode,'Band-pass'),
       ffilter1=shape(lft1,centerfrequency+1,filterwidth,filtershape);
       ffilter2=shape(lft2,length(fy)-centerfrequency+1,filterwidth,filtershape);
       ffilter=[ffilter1,ffilter2];
end
if strcmp(mode,'High-pass')
       centerfrequency=length(x)/2;
       ffilter1=shape(lft1,centerfrequency+1,filterwidth,filtershape);
       ffilter2=shape(lft2,length(fy)-centerfrequency+1,filterwidth,filtershape);
       ffilter=[ffilter1,ffilter2];
end
if strcmp(mode,'Low-pass')
       centerfrequency=0;
       ffilter1=shape(lft1,centerfrequency+1,filterwidth,filtershape);
       ffilter2=shape(lft2,length(fy)-centerfrequency+1,filterwidth,filtershape);
       ffilter=[ffilter1,ffilter2];
end
if strcmp(mode,'Band-reject (notch)')
       ffilter1=shape(lft1,centerfrequency+1,filterwidth,filtershape);
       ffilter2=shape(lft2,length(fy)-centerfrequency+1,filterwidth,filtershape);
       ffilter=1-[ffilter1,ffilter2]; 
end


if length(fy)>length(ffilter), ffilter=[ffilter ffilter(1)];end

ffy=fy.*ffilter;  % Multiply filter by Fourier Transfork of signal
ry=real(ifft(ffy));

if display
   figure, plot(y,'Color','b'); hold on; plot(ry,'Color','r'); 
    
end

py=fy .* conj(fy); % Compute power spectrum
    plotrange=2:length(fy)/2;

   f=((plotrange-1)./range(x));
   
   py=real(py(plotrange));
   
if display>1
   figure
   
   
   
    loglog(f,py,f,max(py).*ffilter(plotrange),'r')

title('POWER SPECTRA:  BLUE = Input signal    RED = Filter')
xlabel('Frequency');
ylabel('Power spectrum (A.U^2/Hz)');
    
end



function g = shape(x,pos,wid,n)
%  shape(x,pos,wid,n) = peak centered on x=pos, half-width=wid
%  x may be scalar, vector, or matrix, pos and wid both scalar
%  Shape is Lorentzian (1/x^2) when n=0, Gaussian (exp(-x^2))
%  when n=1, and becomes more rectangular as n increases.
%  Example: shape([1 2 3],1,2,1) gives result [1.0000    0.5000    0.0625]
if n==0
    g=ones(size(x))./(1+((x-pos)./(0.5.*wid)).^2);
else
    g = exp(-((x-pos)./(0.6.*wid)) .^(2*round(n)));
end

function r = range(arr)
r = max(arr) - min(arr);
