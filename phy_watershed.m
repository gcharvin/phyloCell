function L = phy_watershed(varargin)
%WATERSHED Watershed transform.
%   L = WATERSHED(A) computes a label matrix identifying the watershed
%   regions of the input matrix A.  A can have any dimension.  The elements
%   of L are integer values greater than or equal to 0.  The elements
%   labeled 0 do not belong to a unique watershed region.  These are called
%   "watershed pixels."  The elements labeled 1 belong to the first
%   watershed region, the elements labeled 2 belong to the second watershed
%   region, and so on.
%
%   By default, WATERSHED uses 8-connected neighborhoods for 2-D inputs
%   and 26-connected neighborhoods for 3-D inputs.  For higher
%   dimensions, WATERSHED uses the connectivity given by
%   CONNDEF(NDIMS(A),'maximal').
%
%   L = WATERSHED(A,CONN) computes the watershed transform using the
%   specified connectivity.  CONN may have the following scalar values:
%
%       4     two-dimensional four-connected neighborhood
%       8     two-dimensional eight-connected neighborhood
%       6     three-dimensional six-connected neighborhood
%       18    three-dimensional 18-connected neighborhood
%       26    three-dimensional 26-connected neighborhood
%
%   Connectivity may be defined in a more general way for any dimension by
%   using for CONN a 3-by-3-by- ... -by-3 matrix of 0s and 1s.  The 1-valued
%   elements define neighborhood locations relative to the center element of
%   CONN.  If specified this way, CONN must be symmetric about its center.
%
%   Note
%   ----
%   The watershed transform algorithm used by this function changed in
%   version 5.4 (R2007a) of the Image Processing Toolbox.  The previous 
%   algorithm occasionally produced labeled watershed basins that were not
%   contiguous.  If you need to obtain the same results as the previous
%   algorithm, use the function WATERSHED_OLD.
%
%   Class Support
%   -------------
%   A can be a numeric or logical array of any dimension, and it must be
%   nonsparse.  The output array L is double. 
%
%   Example (2-D)
%   -------------
%   1. Make a binary image containing two overlapping circular objects.
%
%       center1 = -10;
%       center2 = -center1;
%       dist = sqrt(2*(2*center1)^2);
%       radius = dist/2 * 1.4;
%       lims = [floor(center1-1.2*radius) ceil(center2+1.2*radius)];
%       [x,y] = meshgrid(lims(1):lims(2));
%       bw1 = sqrt((x-center1).^2 + (y-center1).^2) <= radius;
%       bw2 = sqrt((x-center2).^2 + (y-center2).^2) <= radius;
%       bw = bw1 | bw2;
%       figure, imshow(bw,'InitialMagnification','fit'), title('bw')
%
%   2. Compute the distance transform of the complement of the binary
%      image. 
%
%       D = bwdist(~bw);
%       figure, imshow(D,[],'InitialMagnification','fit')
%       title('Distance transform of ~bw')
%
%   3. Complement the distance transform, and force pixels that don't
%      belong to the objects to be at -Inf.
%
%       D = -D;
%       D(~bw) = -Inf;
%
%   4. Compute the watershed transform, and display the resulting label
%      matrix as an RGB image.
%
%       L = watershed(D); 
%       rgb = label2rgb(L,'jet',[.5 .5 .5]);
%       figure, imshow(rgb,'InitialMagnification','fit')
%       title('Watershed transform of D')
%
%   Example (3-D)
%   -------------
%   1. Make a 3-D binary image containing two overlapping spheres.
%
%       center1 = -10;
%       center2 = -center1;
%       dist = sqrt(3*(2*center1)^2);
%       radius = dist/2 * 1.4;
%       lims = [floor(center1-1.2*radius) ceil(center2+1.2*radius)];
%       [x,y,z] = meshgrid(lims(1):lims(2));
%       bw1 = sqrt((x-center1).^2 + (y-center1).^2 + ...
%           (z-center1).^2) <= radius;
%       bw2 = sqrt((x-center2).^2 + (y-center2).^2 + ...
%           (z-center2).^2) <= radius;
%       bw = bw1 | bw2;
%       figure, isosurface(x,y,z,bw,0.5), axis equal, title('BW')
%       xlabel x, ylabel y, zlabel z
%       xlim(lims), ylim(lims), zlim(lims)
%       view(3), camlight, lighting gouraud
%
%   2. Compute the distance transform.
%
%       D = bwdist(~bw);
%       figure, isosurface(x,y,z,D,radius/2), axis equal
%       title('Isosurface of distance transform')
%       xlabel x, ylabel y, zlabel z
%       xlim(lims), ylim(lims), zlim(lims)
%       view(3), camlight, lighting gouraud
%
%   3. Complement the distance transform, force nonobject pixels to be
%      -Inf, and then compute the watershed transform.
%
%       D = -D;
%       D(~bw) = -Inf;
%       L = watershed(D);
%       figure, isosurface(x,y,z,L==2,0.5), axis equal
%       title('Segmented object')
%       xlabel x, ylabel y, zlabel z
%       xlim(lims), ylim(lims), zlim(lims)
%       view(3), camlight, lighting gouraud
%       figure, isosurface(x,y,z,L==3,0.5), axis equal
%       title('Segmented object')
%       xlabel x, ylabel y, zlabel z
%       xlim(lims), ylim(lims), zlim(lims)
%       view(3), camlight, lighting gouraud
%
%   See also BWLABEL, BWLABELN, REGIONPROPS, WATERSHED_OLD.

%   Copyright 1993-2006 The MathWorks, Inc.
%   $Revision: 1.10.4.7 $  $Date: 2006/10/14 12:22:54 $

% Input-output specs
% ==================
% A    - full, real, numeric, logical
%        +/- Inf OK, but NaNs not allowed
%        empty OK
%        required
%
% CONN - connectivity; see connectivity spec
%        optional; if not specified = ones(repmat(3,1,ndims(A)))
%
% L      full, double array, same size as A

[A,conn] = parse_inputs(varargin{:});
%M = local_bwlabel(imregionalmin(A,conn),conn);
pix=find(A==min(min(A)));
B=zeros(size(A));
B(pix)=1;
%figure; imshow(B);
M = local_bwlabel(B,conn);
L = watershed_meyer(A,conn,M);

%====================================================================
% The function bwlabel is faster than bwlabeln, but bwlabel only
% supports 2-D inputs, and a conn specified as either 4 or 8.
% local_bwlabel uses bwlabel when possible; otherwise it uses
% bwlabeln.
function L = local_bwlabel(BW, conn)
if ndims(BW) == 2
    if is8Conn(conn)
        labeler = @(BW) bwlabel(BW, 8);
        
    elseif is4Conn(conn)
        labeler = @(BW) bwlabel(BW, 4);
        
    else
        labeler = @(BW) bwlabeln(BW, conn);
    end
else
    labeler = @(BW) bwlabeln(BW, conn);
end

L = labeler(BW);
%--------------------------------------------------------------------

%====================================================================
% is4Conn
% True if the connectivity specifier is either the value 4, or is
% [0 1 0; 1 1 1; 0 1 0].
%
function tf = is4Conn(conn)

tf = isequal(conn, 4) || isequal(logical(conn), logical([0 1 0; 1 1 1; 0 1 0]));
%--------------------------------------------------------------------

%====================================================================
% is8Conn
% True if the connectivity specifier is either the value 8, or is
% ones(3,3).
%
function tf = is8Conn(conn)

tf = isequal(conn, 8) || isequal(logical(conn), true(3,3));
%--------------------------------------------------------------------

%====================================================================
function [A,conn] = parse_inputs(varargin)
  
iptchecknargin(1,2,nargin,mfilename);

A = varargin{1};
iptcheckinput(A,{'numeric' 'logical'}, {'real' 'nonsparse'}, ...
              mfilename, 'A', 1);

if nargin < 2
    conn = conndef(ndims(A), 'maximal');
else
    conn = varargin{2};
    if isa(conn,'strel')
        conn = getnhood(conn);
    end
    iptcheckconn(conn, mfilename, 'CONN', 2);
end
%--------------------------------------------------------------------
