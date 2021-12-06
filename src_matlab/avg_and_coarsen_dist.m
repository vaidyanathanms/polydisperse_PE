%% To plot average distribution and to coarsen the same
% Run out_mwdist.m before using this. Requires file of the form:
% ./../../outfiles/overall/out_mwdist_n_%d_pdi_%g_%s_rcut_%s
% change file name if required.
% 5th column has to be the MW of unique chain list.
% 5th column has to be normalized adsorption probability. See example of an
% input files for details.

clc;
clear;
close all;
format long;

%% Color codes

green = [0 0.5 0.0]; gold = [0.9 0.75 0]; orange = [0.91 0.41 0.17]; brown=[0.6 0.2 0];
pclr = {'m','r',gold,green,orange,'k','b',brown};
lsty = {'-','--',':'};
msty = {'d','s','o','x'};

%% Flags
avg_flag = 1;
plt_flag = 1;

%% Input data
nch_freearr = [32;64;128;150];
casearr  = [1;2;3;4];
pdi_freearr = [1.5];
arch_arr  = {'bl_bl'};
pdigraft  = 1.0;
nfreemons = 30;
ngraft_ch = 32; % Number of graft chains
cutoff = '1.50';
lz = 120;
area = 35^2;
kT = 1;
xlab_req = 1; % Plot with x labels if value is 1
%% Graft details
ncharge_mons_graft = 30; % per graft chain details
ntail_mons_graft = 5;    % per graft chain details
ntot_mons_graft  = ncharge_mons_graft + ntail_mons_graft;
nch_graft = 32;

%% Pre-calculations
rhofree = nch_freearr*nfreemons/(lz*area);
pdigraft_str = num2str(pdigraft,'%1.1f');
num_cases = length(casearr);
max_mw_free = 10*nfreemons; % An approximate max. Will throw error from extract_adschain.m if it is more than this value.

%% For spline fit

pvalarr = [5.8023914813527593E-4;3.5018844316573883E-5;1.1096443546639338E-4;2.47616468322602E-5];
%% Compute average_distribution

if avg_flag
    
    for ncnt = 1:length(nch_freearr) % begin nfree loop
        nval = nch_freearr(ncnt);
        
        for pdi_cntr = 1:length(pdi_freearr) % begin pdi free loop
            ref_pdifree     = pdi_freearr(pdi_cntr);
            pdifree_str     = num2str(ref_pdifree,'%1.1f');
            
            for arch_cnt = 1:length(arch_arr)  % begin arch loop
                dirstr = arch_arr{arch_cnt};
                
                mw_data_arr = zeros(50,3); %value of 50 is approximate. Can weed off zero at the end.
                mw_ref_cntr = 0;
                
                dirname = sprintf('./../../outfiles/overall');
                if ~exist(dirname,'dir')
                    fprintf('%s does not exist\n',dirname);
                    continue
                end
                
                finp_fyle = sprintf('./../../outfiles/overall/out_mwdist_n_%d_pdi_%g_%s_rcut_%s.dat',...
                    nval,ref_pdifree,dirstr,cutoff);
                finp_data = fopen(finp_fyle,'r');
                if finp_data <= 0
                    fprintf('ERROR: %s not found\n', finp_data);
                    continue;
                end
                
                s1 = create_output_dirs('./../../distribution_dir/avg_values');

                fprintf('Analyzing %s\n', finp_fyle);
                fout_data = fopen(sprintf('./../../distribution_dir/avg_values/avg_mwdist_n_%d_pdi_%g_%s_rcut_%s.dat',...
                    nval,ref_pdifree,dirstr,cutoff),'w');
                fprintf(fout_data,'%s\t%s\t%s\t%s\n','MW','Unnorm_probability','Tot_occurences','Norm_probability');
                
                err_flag = 0;
                % Parse and analyze the file.
                while ~feof(finp_data) && err_flag == 0
                    
                    % Start reading file and find header keyword
                    tline = fgetl(finp_data); % get header
                    if ~ischar(tline) || isempty(tline)
                        fprintf('ERROR: Unable to read file %s\n', tline)
                        return;
                    end
                    spl_tline = strtrim(strsplit(strtrim(tline)));
                    find_keyword = -1; % to find "freechainMW" keyword
                    
                    for wordcnt = 1:length(spl_tline)
                        if strcmp(strtrim(spl_tline{wordcnt}),'num_unique_MW')
                            find_keyword = 1; column_num = wordcnt;
                            num_unique_MWs = str2double(strtrim(spl_tline{column_num+1}));
                            clear column_num
                            break;
                        end
                    end
                    
                    % check if the first and fifth column are the MW and
                    % norm_adsorption_prob respectively.
                    tline = fgetl(finp_data); % get column numbers of adsorption prob for each case
                    spl_tline = strtrim(strsplit(strtrim(tline)));
                    if ~strcmp(strtrim(spl_tline{1}),'MW') || ~strcmp(strtrim(spl_tline{5}),'Norm_adsorption_prob')
                        fprintf('ERROR: 1st and 5th column needs to be MW and norm_adsorption_prob respectively: %s\t%s\n',strtrim(spl_tline{1}),strtrim(spl_tline{5}));
                        continue;
                    end
                    
                    
                    for linecnt = 1:num_unique_MWs % Read each case and process
                        tline = fgetl(finp_data);
                        spl_tline = strtrim(strsplit(strtrim(tline)));
                        MW_val = str2double(strtrim(spl_tline{1}));
                        norm_adsorb_val = str2double(strtrim(spl_tline{5}));
                        
                        if MW_val <=0
                            fprintf('ERROR: Unknown mol. wt: %d\n',MW_val);
                            continue;
                        end

                        %find MW_val is already present in the first column mw_data_arr
                        check_flag = ismember(mw_data_arr(:,1),MW_val);
                        if max(check_flag(:,1)) == 0
                            mw_ref_cntr = mw_ref_cntr + 1;
                            mw_data_arr(mw_ref_cntr,1) = MW_val;
                            mw_data_arr(mw_ref_cntr,2) = norm_adsorb_val;
                            mw_data_arr(mw_ref_cntr,3) = 1; %change flag
                        else % if it is already present, increment counter in the third column and the normalized value to second column so that it can be divided at the end
                            index_val = find(mw_data_arr(:,1)==MW_val);
                            if length(index_val) > 1
                                fprintf('ERROR: Multiple occurences of the same MW (%d) found in the consolidated array \n', MW_val);
                                fprintf('Mol wt. array data:\n');
                                fprintf('%d\n',mw_data_arr(:,1));
                                err_flag = 1;
                                break;
                            end
                            
                            if mw_data_arr(index_val,1) ~= MW_val
                                fprintf('ERROR: Not adding the corresponding MWs: %d\t%d\n', mw_data_arr(index_val,1), MW_val);
                                err_flag = 1;
                                break;
                            end
                            
                            mw_data_arr(index_val,2) = mw_data_arr(index_val,2) + norm_adsorb_val;
                            mw_data_arr(index_val,3) = mw_data_arr(index_val,3) + 1;
                        end % end processing one line
                        
                    end % end processing each case
                    
                end % end reading the file for a given architecture (end of while loop)
                
                if err_flag ~= 1
                          
                    for write_cntr = 1:length(mw_data_arr(:,1))
                        
                        if mw_data_arr(write_cntr,3) == 0
                            fprintf('WARNING: Number of elements is less than 50 for this case, RECHCEK \n');
                            continue;
                        end
                        
                        norm_vals = mw_data_arr(write_cntr,2)/mw_data_arr(write_cntr,3);
                        fprintf(fout_data,'%d\t%d\t%d\t%d\n',...
                            mw_data_arr(write_cntr,1),mw_data_arr(write_cntr,2),mw_data_arr(write_cntr,3),norm_vals);
                        
                    end
                    
                end
                
            end % end arch loop
            
        end % end pdi loop
        
    end % end nval loop
    
end % end avg_flag


%% Plot distribution

if plt_flag
    
    % Plot in one figure
    for pdi_cntr = 1:length(pdi_freearr) % begin pdi free loop
        ref_pdifree     = pdi_freearr(pdi_cntr);
        pdifree_str     = num2str(ref_pdifree,'%1.1f');
        
        for arch_cnt = 1:length(arch_arr)  % begin arch loop
            dirstr = arch_arr{arch_cnt};
            
            harch = figure;
            hold on
            box on
            set(gca,'FontSize',16)
            xlabel('$MW$','FontSize',20,'Interpreter','Latex')
            ylabel('$p_{\rm{ads}}(MW)$','FontSize',20,'Interpreter','Latex')
            
            
            for ncnt = 1:length(nch_freearr) % begin nfree loop
                nval = nch_freearr(ncnt);
                
                fylename = sprintf('./../../distribution_dir/avg_values/avg_mwdist_n_%d_pdi_%g_%s_rcut_%s.dat',...
                    nval,ref_pdifree,dirstr,cutoff);
                if exist(fylename,'file') ~=2
                    fprintf('%s does not exist', fylename);
                    continue;
                end
                
                plt_data = importdata(fylename);
                plot(plt_data.data(:,1), plt_data.data(:,4), 'Color',pclr{ncnt},'LineStyle','None',...
                    'Marker',msty{pdi_cntr},'MarkerFaceColor',pclr{ncnt},'MarkerEdgeColor',pclr{ncnt},'MarkerSize',8)
                legendinfo{ncnt} = ['$N_{pa}/N_{pc}$ = ' num2str(nval/ngraft_ch,'%1.1f')];
                
            end % end nfree loop
            ylim([0 1.1]);
            legend(legendinfo,'Interpreter','Latex','FontSize',14,'Location','SouthEast')
            legend boxon
            saveas(harch,sprintf('./../../Figs_paper/avg_dist_%s_%s.png',dirstr,pdifree_str));
            
        end % end arch loop
        
    end % end pdi loop
    
    
    % Plot as subplots
    for pdi_cntr = 1:length(pdi_freearr) % begin pdi free loop
        ref_pdifree     = pdi_freearr(pdi_cntr);
        pdifree_str     = num2str(ref_pdifree,'%1.1f');
        
        for arch_cnt = 1:length(arch_arr)  % begin arch loop
            dirstr = arch_arr{arch_cnt};
            
            harch_subplot = figure;
            hold on
            box on
            axis square
            set(gcf, 'Units', 'normalized', 'Position', [0.5,0.5,0.3,0.1] ) ;
            hs = zeros(4);
            x_min_max_lims = zeros(4,2); %for setting x limits of plot
            y_min_max_lims = zeros(4,2); %for setting y limits of plot
            
            for ncnt = 1:length(nch_freearr) % begin nfree loop
                nval = nch_freearr(ncnt);
                
                fylename = sprintf('./../../distribution_dir/avg_values/avg_mwdist_n_%d_pdi_%g_%s_rcut_%s.dat',...
                    nval,ref_pdifree,dirstr,cutoff);
                if exist(fylename,'file') ~=2
                    fprintf('%s does not exist', fylename);
                    continue;
                end
                
                plt_data = importdata(fylename);
                xdata = plt_data.data(:,1); ydata = plt_data.data(:,4) ;
                hs(ncnt) = subplot(1,4,ncnt);
                plot(plt_data.data(:,1), plt_data.data(:,4), 'Color',pclr{ncnt},'LineStyle','None',...
                    'Marker',msty{pdi_cntr},'MarkerFaceColor',pclr{ncnt},'MarkerEdgeColor',pclr{ncnt},'MarkerSize',8)
                init_guess = [0,0];
                [hval,eval] = fit_to_theory(xdata,ydata,1/kT,init_guess);
                hold on
                %                 F = @(x,xdata)tanh(2*(xdata-x(1))/x(2));
                %                 x0 = [ntot_mons_graft ntot_mons_graft];
                %                 [x,resnorm,~,exitflag,output] = lsqcurvefit(F,x0,plt_data.data(:,1),plt_data.data(:,4));
                %                 plot(0:max(plt_data.data(:,1)),F(x,0:max(plt_data.data(:,1))))
                
                %                 xx = 0:1:1.2*max(plt_data.data(:,1));
                %                 yy = spline(plt_data.data(:,1),plt_data.data(:,4),xx);
                %                 plot(xx,yy,'Color',pclr{ncnt},'LineStyle','-','LineWidth',2)
                
                
                
                % return an array with y~=0 and x<~6 STD of number of chains (surrogate for MW) to draw guideline to eye
                newarrcnt = 1;
                for cpycnt = 1:length(plt_data.data(:,1))
                    if plt_data.data(cpycnt,4) ~= 0 && plt_data.data(cpycnt,1)/nval < 3% very crude way to get the pvalue
                        xalldata(newarrcnt,1) = plt_data.data(cpycnt,1);
                        yalldata(newarrcnt,1) = plt_data.data(cpycnt,4);
                        newarrcnt = newarrcnt+1;
                    end
                end
                
                % Ref:https://www.mathworks.com/help/curvefit/cubic-smoothing-splines.html
                %epsilon = ((xinp(end)-xinp(1))/(numel(xinp)-1))^3/16;
                pval = pvalarr(ncnt); %https://www.mathworks.com/help/curvefit/cubic-smoothing-splines.html
                xxi = (0:0.9*max(xalldata));
                ys = csaps(xalldata,yalldata,pval,xxi,yalldata);
%                 plot(xxi,ys,'Color',pclr{ncnt},'LineStyle','--','LineWidth',2)
                xval = linspace(min(xdata),max(xdata));
                fitdata = (hval*exp(eval*xval))./(1+hval*exp(eval*xval));
%                plot(xval,fitdata,'Color',pclr{ncnt},'LineStyle','--','LineWidth',2)
                if ncnt == 4
                    x3d = xalldata; y3d = yalldata;
                end
                clear xalldata yalldata
                
                x_min_max_lims(ncnt,1) = min(plt_data.data(:,1));
                x_min_max_lims(ncnt,2) = max(plt_data.data(:,1));
                
                y_min_max_lims(ncnt,1) = min(plt_data.data(:,4));
                y_min_max_lims(ncnt,2) = max(plt_data.data(:,4));
                
                lgd = legend(['$n_{pa}/n_{pc}$ = ' num2str(nval/ngraft_ch,'%1.1f')]);
                lgd.FontSize = 25;
                lgd.Interpreter = 'Latex';
                lgd.Location='NorthWest';
                legend boxon
                
            end % end nfree loop
            
            p1 = get(hs(1),'Position'); %[x,y,w,h]
            p2 = get(hs(2),'Position');
            p3 = get(hs(3),'Position');
            p4 = get(hs(4),'Position');
            
            if xlab_req
                set(hs(4),'YTick',[],'XTick',[0 60 120]);
                set(hs(3),'YTick',[],'XTick',[0 60 120]);
                set(hs(2),'YTick',[],'XTick',[0 60 120]);
                set(hs(1),'XTick',[0 60 120]);
            else
                set(hs(4),'YTick',[],'XTick',[0 60 120],'XTickLabel', []);
                set(hs(3),'YTick',[],'XTick',[0 60 120],'XTickLabel', []);
                set(hs(2),'YTick',[],'XTick',[0 60 120],'XTickLabel', []);
                set(hs(1),'XTick',[0 60 120],'XTickLabel', []);
            end
                
            ylim(hs(1),[y_min_max_lims(4,1) y_min_max_lims(4,2)+0.15]);
            ylim(hs(2),[y_min_max_lims(4,1) y_min_max_lims(4,2)+0.15]);
            ylim(hs(3),[y_min_max_lims(4,1) y_min_max_lims(4,2)+0.15]);
            ylim(hs(4),[y_min_max_lims(4,1) y_min_max_lims(4,2)+0.15]);
            
            xmax = 150; % x_min_max_lims(4,2)+0.2*x_min_max_lims(4,2);
            xlim(hs(1),[x_min_max_lims(4,1)-6 xmax]);
            xlim(hs(2),[x_min_max_lims(4,1)-6 xmax]);
            xlim(hs(3),[x_min_max_lims(4,1)-6 xmax]);
            xlim(hs(4),[x_min_max_lims(4,1)-6 xmax]);
            
            p1(3) = 0.775/4;
            p2(3) = 0.775/4;
            p3(3) = 0.775/4;
            p4(3) = 0.775/4;
            
            p1(4) = 0.775;
            p2(4) = 0.775;
            p3(4) = 0.775;
            p4(4) = 0.775;
            
            p1(1) = 0.15;
            p2(1) = p1(1) + p1(3);
            p3(1) = p2(1) + p2(3);
            p4(1) = p3(1) + p3(3);
            
            set(hs(1),'pos', p1,'FontSize',32);
            set(hs(2),'pos', p2,'FontSize',32);
            set(hs(3),'pos', p3,'FontSize',32);
            set(hs(4),'pos', p4,'FontSize',32);
            
            han=axes(harch_subplot,'visible','off');
            han.XLabel.Visible='on';
            han.YLabel.Visible='on';
            ylabel(hs(1),'$p_{\rm{ads}}$','FontSize',40,'Interpreter','Latex');
            if xlab_req
                xlb = xlabel(hs(2),'$N$','FontSize',40,'Interpreter','Latex');
            end
            
            saveas(harch_subplot,sprintf('./../../Figs_paper/avg_dist_%s_%s.png',dirstr,pdifree_str));
            
        end % end arch loop
        
    end % end pdi loop
    
    
end % end plt_flag loop
