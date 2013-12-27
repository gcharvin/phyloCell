% change from cells in struct format to cells in class format
% from old type to new type
function cells_class=phy_struct2class(cells_struct)
nx=32;
cells_class=phy_Object;
if isfield(cells_struct,'n')
    for i=1:cells_struct.n
        n=i-1;
        x=cells_struct.x(n*(nx+1)+1:(n+1)*(nx+1));
        y=cells_struct.y(n*(nx+1)+1:(n+1)*(nx+1));
        cells_class(i).x=x;
        cells_class(i).y=y;
        cells_class(i).ox=cells_struct.ox(i);
        cells_class(i).oy=cells_struct.oy(i);
        cells_class(i).area=cells_struct.area(i);
        if isfield(cells_struct,'displayMap') && length(cells_struct.displayMap)==cells_struct.n % if already named
            cells_class(i).n=cells_struct.displayMap(i); %kep the old name
        else
            cells_class(i).n=i;
        end
    end
end
