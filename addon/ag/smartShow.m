function smartShow(im)
    persistent xlims;
    persistent ylims;
    
    if ~isempty(get(gcf, 'CurrentAxes'))
        xlims = xlim;
        ylims = ylim;
    end
    
    imshow(im);
    
    if ~isempty(xlims)
        xlim(xlims)
        ylim(ylims)
    end
end