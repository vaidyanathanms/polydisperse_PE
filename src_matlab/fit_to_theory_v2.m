function vaidya_theory_v2
clc
close all
set(0,'defaulttextinterpreter','latex')

global xmin xmax
xmin = 0; %smallest N value for plotting purposes
xmax = 160; % biggest N value for plotting purposes

fprintf('======================================\n')
fprintf('Analyzing the block-block data \n')
fprintf('======================================\n')

analyze_data('blbl')


fprintf('======================================\n')
fprintf('Analyzing the alt-alt data \n')
fprintf('======================================\n')

analyze_data('alal')



function analyze_data(sim_type)
global xmin xmax

%Vaidya colors
green = [0 0.5 0.0]; 
gold = [0.9 0.75 0];
pclr = {'m','r',gold,green};

%upper bounds for doing the regression
Nmax = [25,40,80,80];


%read in the data from the Figure 7 file
fprintf('Reading in the data file \n')
[z47_full,z40_full,z20_full,z10_full] = get_data_from_fig(sim_type);

%trim out the zero and one data that will give singularities
%ncut tells you how many data points were removed. 
%the lower and upper bounds are where to cut the data so that
%we are only looking at the bulk of the distribution
p_lower = 0.0001;
p_upper = 0.9999;
fprintf('Removing data points with p = 0 or p = 1 \n')
[z47,ncut] = trim_zeros(z47_full,p_lower,p_upper);
fprintf('\t Removed %2d points for n = 4.7\n',ncut)
[z40,ncut] = trim_zeros(z40_full,p_lower,p_upper);
fprintf('\t Removed %2d points for n = 4.0\n',ncut)
[z20,ncut] = trim_zeros(z20_full,p_lower,p_upper);
fprintf('\t Removed %2d points for n = 2.0\n',ncut)
[z10,ncut] = trim_zeros(z10_full,p_lower,p_upper);
fprintf('\t Removed %2d points for n = 1.0\n',ncut)

%make a plot of all data with regression
fprintf('Making a plot of all of the data \n')
plot_semilog(sim_type,z10,z20,z40,z47,pclr,Nmax)

%plot each of the data sets on their own w/o log formatting
fprintf('Making a plot of each data set \n')
plot_single_data(sim_type,z10_full,z10,pclr,Nmax,1)
plot_single_data(sim_type,z20_full,z20,pclr,Nmax,2)
plot_single_data(sim_type,z40_full,z40,pclr,Nmax,3)
plot_single_data(sim_type,z47_full,z47,pclr,Nmax,4)

function plot_single_data(sim_type,zfull,z,pclr,Nmax,panel_num)
%plot a single case with a distribution
%z = trimmed data
%zfull = the full data set, not trimmed
%panel_num = the color choice/Nmax for that panel
global xmax xmin
g = figure;
hold on 
plot_data(zfull,pclr(panel_num),'d');
[a,b] = regress_theory(z,Nmax(panel_num));
xreg = linspace(1,xmax,101);
w = a*exp(b*xreg);
yreg = w./(1+w);
plot(xreg,yreg,'-','Color',pclr{panel_num},'LineWidth',2)
hold off
xlabel_str = '$N$';
ylabel_str = '$p_{\rm ads}$';
xlims = [0,120];
ylims = [0,1];
xscale = 'linear';
yscale = 'linear';
format_plot(xlabel_str,ylabel_str,xlims,ylims,xscale,yscale)
ytick=0:0.2:1;
YTickLabels = cellstr(num2str(ytick(:)));
set(gca,'YTick',ytick,'YTickLabel',YTickLabels,'TickLabelInterpreter','tex')
filename = strcat(sim_type,'_panel',num2str(panel_num),'_pads.eps');
saveas(g,strcat('../../Figs_paper/',filename),'epsc')


function plot_semilog(sim_type,z10,z20,z40,z47,pclr,Nmax)
global xmin xmax
%plot all the data together and output regression coefficients

%set the bounds for N values for plotting here so that we can find them
ymin = 1e-3; %had to adjust these numbers to get a decent value
ymax = 200; %had to adjust this number to get a good fit

%This is just a vector for my plot formatting command
xlims = [xmin,xmax]; 
ylims = [ymin,ymax];



g = figure;
hold on 
plot_data(z10,pclr(1),'d');
plot_data(z20,pclr(2),'d');
plot_data(z40,pclr(3),'d');
plot_data(z47,pclr(4),'d');
%do the linear regression and write the results to file
fileout = strcat('../../outfiles/overall/',sim_type,'_regression.dat');
fid = fopen(fileout,'w');
fprintf(fid,strcat('Regression data for \t ',sim_type,'\n'));
fprintf(fid,'Case \t Nmax \t a \t\t b \n');
    [a,b] = regress_theory(z10,Nmax(1));
    fprintf(fid,'z10 \t %2d \t %6.4f \t %6.4f \n',Nmax(1),a,b);
    plot_regression(a,b,pclr(1))
    [a,b] = regress_theory(z20,Nmax(2));
    fprintf(fid,'z20 \t %2d \t %6.4f \t %6.4f \n',Nmax(2),a,b);
    plot_regression(a,b,pclr(2))
    [a,b] = regress_theory(z40,Nmax(3));
    fprintf(fid,'z40 \t %2d \t %6.4f \t %6.4f \n',Nmax(3),a,b);
    plot_regression(a,b,pclr(3))
    [a,b] = regress_theory(z47,Nmax(4));
    fprintf(fid,'z10 \t %2d \t %6.4f \t %6.4f \n',Nmax(4),a,b);
    plot_regression(a,b,pclr(4))
    fclose(fid);
hold off
L10 = '$n_{pa}/n_{pc} = 1.0$';
L20 = '$n_{pa}/n_{pc} = 2.0$';
L40 = '$n_{pa}/n_{pc} = 4.0$';
L47 = '$n_{pa}/n_{pc} = 4.7$';
legend({L10,L20,L40,L47},'Interpreter','Latex',...
    'Location','SouthEast','FontSize',18)
xlabel_str = '$N$';
ylabel_str = '$p_{\rm ads}/(1-p_{\rm ads})$';
xscale = 'linear';
yscale = 'log';
format_plot(xlabel_str,ylabel_str,xlims,ylims,xscale,yscale)
ytick=10.^(-2:2:2);
YTickLabels = cellstr(num2str(round(log10(ytick(:))), '10^{%d}'));
set(gca,'YTick',ytick,'YTickLabel',YTickLabels,'TickLabelInterpreter','tex')
filename = strcat('../../Figs_Paper/',sim_type,'_all_data.eps');
saveas(g,filename,'epsc')


function [a,b] = regress_theory(z,Nmax)
%extract the data for N < Nmax
k = 0;
for i = 1:length(z)
    if z(1,i) <= Nmax
        k = k+1;
        x(k) = z(1,i);
        p(k) = z(2,i);
    end
end

%do linear regression of log p/(1-p) versus N for the data set
y = p./(1-p); %the function from Dave's theory
log_y = log(y); %log of the function from Dave's theory
P = polyfit(x,log_y,1); %log y = P(1)*x + P(2)
a = exp(P(2)); %prefactor of exponential 
b = P(1); %exponential term

function plot_regression(a,b,pclr)
global xmin xmax
%plot the result
xreg = [xmin,xmax]; %bounds for N on regression
yreg = [a*exp(b*xmin), a*exp(b*xmax)]; %bounds converted to linear form
plot(xreg,yreg,'-','Color',pclr{1},'LineWidth',2)


function plot_data(z,pclr,msty)
%plots the function from Dave's theory with the format string
x = z(1,:);
p = z(2,:);
y = p./(1-p);
plot(x,y,'Marker',msty,'MarkerFaceColor',pclr{1},...
    'MarkerEdgeColor',pclr{1},'LineStyle','none','LineWidth',2)



function [ztrim,ncut] = trim_zeros(z,ymin,ymax)
%cut out data that have very low p or p close to 1 because it will cause
%Dave's function to diverge. We only want to fit the data that are in the
%bulk of the distribution. You can change the cutoffs to trim the data to
%be more or less to the center

%remove anything that has a zero or a 1 value
npts = size(z,2);
k = 0;
for i = 1:npts
    if z(2,i) > ymin && z(2,i) < ymax
        k = k + 1;
        ztrim(:,k) = z(:,i);
    end
end
ncut = npts - k;    

function [z47,z40,z20,z10] = get_data_from_fig(sim_type)
%reads the data from Matlab into vectors
%first row is N
%second row is p_ads

filename = strcat('../../Figs_paper/fig7_',sim_type,'_distribution.fig');

%open the figure to get the data
fig = openfig(filename);

%get the xdata
dataObjs = findobj(fig,'-property','Xdata');
x47 = dataObjs(1).XData;
x40 = dataObjs(2).XData;
x20 = dataObjs(3).XData;
x10 = dataObjs(4).XData;

%get the ydata
dataObjs = findobj(fig,'-property','Ydata');
y47 = dataObjs(1).YData;
y40 = dataObjs(2).YData;
y20 = dataObjs(3).YData;
y10 = dataObjs(4).YData;

close %close the figure window

%pack everything into first row = x, second row = y
z47 = [x47; y47];
z40 = [x40; y40];
z20 = [x20; y20];
z10 = [x10; y10];

function format_plot(xlabel_str,ylabel_str,xlims,ylims,xscale,yscale)
%This is just a script I use to make nice figures
    ax = gca;
    ax.FontSize = 18;
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


