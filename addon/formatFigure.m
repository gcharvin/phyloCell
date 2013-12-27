function formatFigure(fontSize,lineSize,color,xlime,ylime,handle)

if nargin<6
    handle=gca;
end

set(handle,'LineWidth',lineSize);
set(handle,'FontSize',fontSize);

%xlabel('Duration of heat-shock (min)');
%ylabel('Apmlitude of HSP104 response (normalized)');

xlab=get(handle,'XLabel');
set(xlab,'FontSize',fontSize);
ylab=get(handle,'YLabel');
set(ylab,'FontSize',fontSize);

xlim(handle,xlime);
ylim(handle,ylime);

if numel(color)~=0
    plots=get(handle,'Children');
    for  k=1:length(plots);
        h=plots(k);
        set(h,'LineWidth',lineSize);
        set(h,'Color',color(k,:));
    end
end

