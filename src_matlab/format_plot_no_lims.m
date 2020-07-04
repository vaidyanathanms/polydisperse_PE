function format_plot_no_lims(xlabel_str,ylabel_str,xscale,yscale)    
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
    %axis labels
    ax.XLabel.String = xlabel_str;
    ax.YLabel.String = ylabel_str;