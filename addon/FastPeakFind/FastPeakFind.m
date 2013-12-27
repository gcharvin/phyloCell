function  cent=FastPeakFind(d, threshold, filt ,edg, fid)
% Analyzes noisy 2D images , finds x-y positions of peaks to 1 pixel accuracy
% The code was meant to be as fast as possible, so I kept it pretty basic.
% The code assumes that the peaks are relatively sparse, test whether there
% is too much pile up and set threshold or user defined filter accordingly.
%
% Inputs:
%   d           The 2D data raw image - assumes a Double\Single-precision
%               floating-point or unit16 array. Please note that the code
%               casts the raw image to uint16, if your image is 8-bit depth,
%               I recommend to change all the places that uint16 is used to
%               uint8 for faster run times. If your image dynamic range is between 0 and 1, I
%               multiplied to fit uint16. This might not be optimal for
%               generic use, so modify according to your needs.
%   threshold   A number between 0 and max(raw_image(:)) to remove  background
%   filt        A filter matrix used to smooth the image. The filter size
%               should correspond the characteristic size of the peaks
%   edg         A number>1 for skipping the first few and the last few 'edge' pixels
%   fid         In case the user would like to save the peak positions to a
%               file, the code assumes a "fid = fopen([filename], 'w+');" line in the script that calls this functions
%
% Output:
%   cent        a 1xN vector of coordinates of peaks (x1,y1,x2,y2,...)
%
%   Example:
%
%   p=FastPeakFind(image);
%   imagesc(image); hold on
%   plot(p(2:2:end),p(1:2:end),'r+')
%
%
%   Nate (nate2718281828@gmail.com)
%   Ver 1.41 , Date: Feb 10th 2013
%% defaults
if (nargin < 1)
    d=uint16(conv2(reshape(single( 2^14*(rand(1,1024*1024)>0.99995) ),[1024 1024]) ,fspecial('gaussian', 15,3),'same')+2^8*rand(1024));
    imagesc(d);
end

if ndims(d)>2 %I added this in case one uses imread (JPG\PNG\...).
    d=uint16(rgb2gray(d));
end

if isfloat(d) %For the case the input image is double, casting to uint16 keeps enough dynamic range while speeds up the code.
    if max(d(:))<=1
        d =  uint16( d.*2^16./(max(d(:))));
    else
        d = uint16(d);
    end
end

if (nargin < 2)
    threshold = (max([min(max(d,[],1))  min(max(d,[],2))])) ;
end

if (nargin < 3)
    filt = (fspecial('gaussian', 7,1)); %if needed modify the filter according to the expected peaks sizes
end

if (nargin < 4)
    edg =3;
end

if (nargin < 5)
    savefileflag = false;
else
    savefileflag = true;
end

%%
if any(d(:))  ; %for the case of non zero raw image

    %d = medfilt2(d,[3,3]);
    
    % apply threshold
    d=d.*uint16(d>threshold);
    
    if any(d(:))   ; %for the case of the image is still non zero
        
        % smooth image
        d=conv2(single(d),filt,'same') ;
        
        % Apply again threshold (and change if needed according to SNR)
        
        d=d.*(d>0.9*threshold);
        
        
        % peak find - using the local maxima approach - 1 pixle resolution
        % d will be noisy on the edges, since no hits are expected there anyway we'll skip 'edge' pixels.
        edg=2;
        [x y]=find(d(edg:size(d,1)-edg,edg:size(d,2)-edg));
        
        % initialize peak find outputs
        cent=[];%
        x=x+edg;
        y=y+edg;
        for j=1:length(y)
            if (d(x(j),y(j))>=d(x(j)-1,y(j)-1 )) &&...
                    (d(x(j),y(j))>d(x(j)-1,y(j))) &&...
                    (d(x(j),y(j))>=d(x(j)-1,y(j)+1)) &&...
                    (d(x(j),y(j))>d(x(j),y(j)-1)) && ...
                    (d(x(j),y(j))>d(x(j),y(j)+1)) && ...
                    (d(x(j),y(j))>=d(x(j)+1,y(j)-1)) && ...
                    (d(x(j),y(j))>d(x(j)+1,y(j))) && ...
                    (d(x(j),y(j))>=d(x(j)+1,y(j)+1));
                
                %All these alternatives were slower...
                %if all(reshape( d(x(j),y(j))>=d(x(j)-1:x(j)+1,y(j)-1:y(j)+1),9,1))
                %if  d(x(j),y(j)) == max(max(d((x(j)-1):(x(j)+1),(y(j)-1):(y(j)+1))))
                %if  d(x(j),y(j))  == max(reshape(d(x(j),y(j))  >=  d(x(j)-1:x(j)+1,y(j)-1:y(j)+1),9,1))
                
                
                %cent(x(j),y(j))=cent(x(j),y(j))+1; % if ones want a matrix output
                cent = [cent ; x(j) ; y(j)];
            end
        end
        
        if savefileflag
            % previous version used dlmwrite, which can be slower than  fprinf
            %             dlmwrite([filename '.txt'],[cent],   '-append', ...
            %                 'roffset', 0,   'delimiter', '\t', 'newline', 'unix');+
            
            fprintf(fid, '%f ', cent(:));
            fprintf(fid, '\n');
            
        end
        
    else % in case image after threshold is all zeros
        cent=[];
        return
    end
    
else % in case raw image is all zeros (dead event)
    cent=[];
    return
end

if (nargin < 1); colormap(bone);hold on; plot(cent(2:2:end),cent(1:2:end),'rs');hold off; end