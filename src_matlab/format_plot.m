function format_plot(xlabel_str,ylabel_str,xlims,ylims,xscale,yscale)    
    ax = gca;
    ax.FontSize = 14;
    ax.Box = 'on';
    ax.YMinorTick = 'on';
    ax.XMinorTick = 'on';
    ax.TickLength = [0.025,0.05];
    ax.XLabel.FontSize = 24;
    ax.YLabel.FontSize = 24;
    ax.XScale = xscale;
    ax.YScale = yscale;
    ax.PlotBoxAspectRatio = [1.618,1,1];
    ax.XColor = 'k';
    ax.YColor = 'k';
    ax.XLim = xlims;
    ax.YLim = ylims;
    %axis labels
    ax.XLabel.String = xlabel_str;
    ax.YLabel.String = ylabel_str;