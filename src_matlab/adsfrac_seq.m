% %% To compute average adsorbed fraction of chains
% 1) Change path to filename (search for keyword filename) and other places
% wherever necessary.

clear;
clc;
close all;
format long;

%% Color codes
green = [0 0.5 0.0]; gold = [0.9 0.75 0]; orange = [0.91 0.41 0.17];brown = [0.2 0 0];
pclr = {'m',brown,green,'k','m', gold};
lsty = {'-','--',':'};
msty = {'d','s','o','x'};

%% Inputs
nfreearr = [150];
casearr  = [1,2,3,4];
pdi_freearr = [1,1.5];
arch_arr = {'bl_bl','bl_al','al_bl','al_al'};
leg_arr  = {'Block-Alter','Block-Alter','Alter-Block','Alter-Alter'}; % ALWAYS CHECK for correspondence with arch_arr
pdigraft = 1.0;
nmonfree = 30; nmongraft = 30; ngraft = 32;
cutoff = '1.50';
lz = 120; area=35^2;
set_tmax = 3e7; % maximum timestep for analysis;
set_tmin = 1e7; % minimum timestep for analysis;

%% Input flags
ttestflag = 1; % write to ttest_dir the individual cases
plotads   = 1; % plot fads as a function of PDI

%% Zero arrays
avg_across_cases = zeros(length(nfreearr),length(arch_arr),length(pdi_freearr));

%% Pre-calculations
rhofree = nfreearr*30/(lz*area);
pdigraft_str = num2str(pdigraft,'%1.1f');

%% Main Analysis

% Create output directories
s1 = create_output_dirs('./../../seq_analysis');
s2 = create_output_dirs('./../../seq_analysis/overall');
s3 = create_output_dirs('./../../ttest_dir');
s4 = create_output_dirs('./../../ttest_dir/overall');

% Create consolidated list
fout_cons = fopen(sprintf('./../../seq_analysis/overall/adsorbed_chain_consolidated_rcut_%s.dat',...
    cutoff),'w');
fprintf(fout_cons,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n','PDI_free','N_f','Arch',...
    'Case_num','min_time','max_time','numsample_pts','avg_fraction');


for pdi_cntr = 1:length(pdi_freearr) % begin pdi free loop
    pdifree     = pdi_freearr(pdi_cntr);
    pdifree_str = num2str(pdifree,'%1.1f');
    
    %zero arrays for averages across cases
    casecntr_arr  = zeros(length(nfreearr),length(arch_arr));
    nadschain_all = zeros(length(nfreearr),length(arch_arr));
    totsamples    = zeros(length(nfreearr),length(arch_arr));
    
    % Create average across all cases
    fout_avg = fopen(sprintf('./../../seq_analysis/overall/adsorbed_chain_ave_allcases_rcut_%s_pdifree_%g.dat',...
        cutoff,pdifree),'w');
    fprintf(fout_avg,'%s\t%s\t%s\t%s\t%s\t%s\n','N_f','Arch','ncases','numsample_pts','avg_fraction','StdErrMean');
    
    % Create case-based avg seq_analysis for SI data
    fout_sidata = fopen(sprintf('./../../seq_analysis/overall/sidata_adsorbed_chain_rcut_%s_pdifree_%g.dat',...
        cutoff,pdifree),'w');
    fprintf(fout_sidata,'%s\t%s\t%s\t%s\t%s\n','\DJ$_{\rm{ideal}}$','$N_{pa}$','Arch',...
        'Case \#','avg_fraction');
    
    
    for ncnt = 1:length(nfreearr) % begin nfree loop
        nval = nfreearr(ncnt);
        
        for arch_cnt = 1:length(arch_arr)  % begin arch loop
            dirstr = arch_arr{arch_cnt};
            
            dirname = sprintf('./../../seq_analysis/outresults_dir_n_%d/%s/pdifree_%s_pdigraft_%s',...
                nval,dirstr,pdifree_str,pdigraft_str);
            
            if ~exist(dirname,'dir')
                fprintf('%s does not exist\n',dirname);
                continue
            end
            
            % Create case-based avg seq_analysis
            s5 = create_output_dirs(sprintf('./../../seq_analysis/n_%d',nval));
            fout_case = fopen(sprintf('./../../seq_analysis/n_%d/adsorbed_chain_rcut_%s_pdifree_%g_%s.dat',...
                nval,cutoff,pdifree,dirstr),'w');
            fprintf(fout_case,'%s\t%s\t%s\t%s\t%s\t%s\n','N_f','Casenum',...
                'Mintime','Maxtime','numavgpoints','avg_fraction');
            
            
            if ttestflag %write all cases into separate folders in ttest_dir
                s5 = create_output_dirs(sprintf('./../../ttest_dir/n_%d',nval));
                fout_ttest = fopen(sprintf('./../../ttest_dir/n_%d/adsfrac_rcut_%s_pdifree_%g_arch_%s.dat',...
                    nval,cutoff,pdifree,dirstr),'w');
                fprintf(fout_ttest,'%s\n','avg_fraction: #of columns correspond to the number of cases');
                
            end
            
            
            avgcase_store = zeros(length(casearr),1); % store each average value
            
            for casecntr = 1:length(casearr) % begin case loop
                casenum = casearr(casecntr);
                
                dirname = sprintf('./../../seq_analysis/outresults_dir_n_%d/%s/pdifree_%s_pdigraft_%s/Case_%d',...
                    nval,dirstr,pdifree_str,pdigraft_str,casenum);
                
                if ~exist(dirname,'dir')
                    fprintf('%s does not exist\n',dirname);
                    continue
                end
                
                fprintf('Analyzing pdi/nfree/arch/case: %g\t%d\t%s\t%d\n', pdifree,nval,dirstr,casenum);
                
                %check if file type exists
                ads_prefix = sprintf('adsfracchain_rcut_%s_config_*.lammpstrj',cutoff);
                ads_fylelist = dir(strcat(dirname,'/',ads_prefix));
                if min(size(ads_fylelist)) == 0
                    fprintf('No files/Empty files are found for %s\n',ads_prefix);
                    continue;
                end
                
                nfyles = numel(ads_fylelist); %number of files of the type
                reordered_ads_list = renumber_files(ads_fylelist,nfyles); % reorder file names to avoid double counting
                
                sum_across_files = 0; tot_cntr_across_files = 0; mintime = 10^10; maxtime = 0;
                mintstep = 0;
                for fylcnt = 1:nfyles % begin running through all files of the given type
                    ads_fylename = strcat(dirname,'/',reordered_ads_list{fylcnt});
                    if exist(ads_fylename,'file') ~= 2
                        fprintf('%s does not exist/empty file\n',ads_fylename);
                        continue;
                    elseif struct(dir(ads_fylename)).bytes == 0
                        fprintf('Empty file: %s \n',ads_fylename);
                        continue;
                    end
                    
                    fprintf('Analyzing %s\n', ads_fylename);
                    data = importdata(ads_fylename);
                    lendata = length(data(:,1));
                    %average adsorption values
                    if min(data(:,1)) > set_tmin
                        for minindcnt = 1:lendata %avoid double counting
                            if data(minindcnt,1) > mintstep
                                minindana = minindcnt; %minimum value at which the trajectories are separate
                                mintstep = max(data(:,1)); %new value will be the maximum value of this file
                                break;
                            end
                        end
                        
                        nads_fracchain = sum(data(minindana:lendata,3));
                        sum_across_files = sum_across_files + nads_fracchain;
                        tot_cntr_across_files = tot_cntr_across_files + length(data(minindana:lendata,3));
                        
                        %find minimum and maximum time
                        if min(data(:,1)) < mintime
                            mintime = min(data(:,1));
                        end
                        if maxtime < max(data(:,1))
                            maxtime = max(data(:,1));
                        end
                        
                    end
                    
                end % end summing adsfrac across all files for a given case
%                 fprintf('%g\t%g\n',sum_across_files,tot_cntr_across_files);
                if tot_cntr_across_files ~= 0
                    avg_for_each_casenum = sum_across_files/tot_cntr_across_files;
                else
                    avg_for_each_casenum = 0;
                end
                
                fprintf(fout_sidata,'%s\t%d\t%s\t%d\t%g\n',pdifree_str,nval,dirstr,...
                    casenum,avg_for_each_casenum);
                
                if maxtime > set_tmax
                    fprintf(fout_case,'%d\t%d\t%d\t%d\t%d\t%g\n',nval,casenum,...
                        mintime,maxtime,tot_cntr_across_files,avg_for_each_casenum);
                    fprintf(fout_cons,'%g\t%d\t%s\t%d\t%d\t%d\t%d\t%g\n',pdifree,nval,...
                        dirstr,casenum,mintime,maxtime,tot_cntr_across_files,avg_for_each_casenum);
                else
                    fprintf(fout_case,'%d\t%d\t%d\t%d\t%d\t%g\t%s\n',nval,casenum,...
                        mintime,maxtime,tot_cntr_across_files,avg_for_each_casenum,'Incomplete sampling');
                    fprintf(fout_cons,'%g\t%d\t%s\t%d\t%d\t%d\t%d\t%g\t%s\n',pdifree,nval,...
                        dirstr,casenum,mintime,maxtime,tot_cntr_across_files,avg_for_each_casenum,'Incomplete sampling');
                end
                
                %save it to overall arrays
                casecntr_arr(ncnt,arch_cnt)  = casecntr_arr(ncnt,arch_cnt) + 1;
                nadschain_all(ncnt,arch_cnt) = nadschain_all(ncnt,arch_cnt) + avg_for_each_casenum;
                totsamples(ncnt,arch_cnt)    = totsamples(ncnt,arch_cnt) + tot_cntr_across_files;
                avgcase_store(casecntr)      = avg_for_each_casenum;
                
                
                if ttestflag % add to corresponding ttest file
                    fprintf(fout_ttest,'%g\t',avg_for_each_casenum);
                end
                
            end % end case loop
            
            fclose(fout_ttest);
            fclose(fout_case);
            
            
            
            if casecntr_arr(ncnt,arch_cnt) > 1
                
                nonzerovals = avgcase_store(avgcase_store~=0);
                stderr_val = std(nonzerovals)/sqrt(length(nonzerovals));
                avg_across_cases(ncnt,arch_cnt,pdi_cntr) = nadschain_all(ncnt,arch_cnt)/length(nonzerovals);
                
            else
                
                stderr_val = -1; % Invalid for just one case
                
            end
            
            fprintf(fout_avg,'%d\t%s\t%d\t%d\t%g\t%g\n',nval,dirstr,casecntr_arr(ncnt,arch_cnt),totsamples(ncnt,arch_cnt),avg_across_cases(ncnt,arch_cnt,pdi_cntr),stderr_val);
            
        end % end arch loop
        
    end % end nfree loop
    
    fclose(fout_avg);
    
    if plotads
        h1 = figure;
        hold on
        box on
        set(gca,'FontSize',16)
        xlabel('$N_{pa}/N_{pc}$','FontSize',20,'Interpreter','Latex')
        ylabel('$f_{ads}$','FontSize',20,'Interpreter','Latex')
        
        for plcnt = 1:length(arch_arr)
            plot(nfreearr/ngraft,avg_across_cases(:,plcnt,pdi_cntr),'color',pclr{plcnt},...
                'Marker',msty{plcnt},'MarkerSize',8,'MarkerFaceColor',pclr{plcnt},'LineStyle','None')
            legendinfo{plcnt} = leg_arr{plcnt};
        end
        %overlay y = x line
        %         x = 0:0.5:5; y = x;
        %         plot(x,y,'LineWidth',2,'Color',orange,'LineStyle','--')
        %         ylim([min(min(avg_across_cases(:,:,pdi_cntr))) 1.2*max(max(avg_across_cases(:,:,pdi_cntr)))])
        legend(legendinfo,'FontSize',16,'Location','Best')
        legend boxoff
        saveas(h1,sprintf('./../../Figs_paper/adsorbchain_rcut_%s_pdi_%g.png',cutoff,pdifree));
        clear legendinfo
    end
    
end % end pdi free loop

fclose(fout_cons);


%% Plot as a function of PDI for each nval

if plotads
    
    for ncnt = 1:length(nfreearr) % begin nfree loop
        nval = nfreearr(ncnt);
        
        h1 = figure;
        hold on
        box on
        set(gca,'FontSize',16)
        xlabel('PDI','FontSize',20,'Interpreter','Latex')
        ylabel('$f_{ads}$','FontSize',20,'Interpreter','Latex')
        
        for plcnt = 1:length(arch_arr)  % begin arch loop
            
            dirstr = arch_arr{plcnt};
            data_to_plot = zeros(length(pdi_freearr),1);
            for pdicntr = 1:length(pdi_freearr)
                data_to_plot(pdicntr) = avg_across_cases(ncnt,plcnt,pdicntr);
            end
            data_to_plot(isnan(data_to_plot)) = 0; %to avoid NaNs
            
            plot(pdi_freearr,data_to_plot,'color',pclr{plcnt},...
                'Marker',msty{plcnt},'MarkerSize',8,'MarkerFaceColor',pclr{plcnt},'LineStyle','None')
            legendinfo{plcnt} = leg_arr{plcnt};
            
        end
        
        legend(legendinfo,'FontSize',16,'Location','Best')
        legend boxoff
        saveas(h1,sprintf('./../../Figs_paper/adsorbchain_rcut_%s_nval_%g.png',cutoff,nval));
        clear legendinfo
        
    end
    
end


