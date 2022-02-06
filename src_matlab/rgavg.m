%% To compute radius of gyration of grafted monomers

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
nfreearr = [16,32,64,128,150];
casearr  = [1,2,3,4];
pdi_freearr = [1,1.5];
arch_arr = {'bl_bl','al_al'};
leg_arr  = {'Block-Block','Alter-Alter'}; % ALWAYS CHECK for correspondence with arch_arr
pdigraft = 1.0;
nmonfree = 30; nmongraft = 30; ngraft = 32;
cutoff = '1.50';
lz = 120; area=35^2;
set_tmax = 2.5e7; % maximum timestep for analysis;
set_tmin = 1e7; % minimum timestep for analysis;

%% Input flags
plotrg   = 1; % plot fads as a function of PDI

%% Zero arrays
avg_across_cases = zeros(length(nfreearr),length(arch_arr),length(pdi_freearr));

%% Pre-calculations
rhofree = nfreearr*30/(lz*area);
pdigraft_str = num2str(pdigraft,'%1.1f');

%% Main Analysis

s1 = create_output_dirs('./../../rganalysis');
s2 = create_output_dirs('./../../rganalysis/overall');

% Create consolidated list
fout_cons = fopen(sprintf('./../../rganalysis/overall/rgavg_consolidated_rcut_%s.dat',...
    cutoff),'w');
fprintf(fout_cons,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n','PDI_free','N_f','Arch',...
    'Case_num','min_time','max_time','numsample_pts','Rgsq_avg');


for pdi_cntr = 1:length(pdi_freearr) % begin pdi free loop
    pdifree     = pdi_freearr(pdi_cntr);
    pdifree_str = num2str(pdifree,'%1.1f');
    
    %zero arrays for averages across cases
    casecntr_arr  = zeros(length(nfreearr),length(arch_arr));
    rg_all = zeros(length(nfreearr),length(arch_arr));
    totsamples    = zeros(length(nfreearr),length(arch_arr));
    
    % Create average across all cases
    fout_avg = fopen(sprintf('./../../rganalysis/overall/rgavg_allcases_rcut_%s_pdifree_%g.dat',...
        cutoff,pdifree),'w');
    fprintf(fout_avg,'%s\t%s\t%s\t%s\t%s\t%s\n','N_f','Arch','ncases','numsample_pts','Rgsq_avg','StdErrMean');
    
    % Create case-based avg outfiles for SI data
    fout_sidata = fopen(sprintf('./../../monads/overall/rgavg_rcut_%s_pdifree_%g.dat',...
        cutoff,pdifree),'w');
    fprintf(fout_sidata,'%s\t%s\t%s\t%s\t%s\n','\DJ$_{\rm{ideal}}$','$N_{pa}$','Arch',...
        'Case \#','Rgsq_avg');
    
    
    for ncnt = 1:length(nfreearr) % begin nfree loop
        nval = nfreearr(ncnt);
        
        for arch_cnt = 1:length(arch_arr)  % begin arch loop
            dirstr = arch_arr{arch_cnt};
            
            dirname = sprintf('./../../rganalysis/outresults_dir_n_%d/%s/pdifree_%s_pdigraft_%s',...
                nval,dirstr,pdifree_str,pdigraft_str);
            
            if ~exist(dirname,'dir')
                %                 fprintf('%s does not exist\n',dirname);
                continue
            end
            
            s3 = create_output_dirs(sprintf('./../../rganalysis/n_%d',nval));
            % Create case-based avg outfiles
            fout_case = fopen(sprintf('./../../rganalysis/n_%d/rgavg_pdifree_%g_%s.dat',...
                nval,pdifree,dirstr),'w');
            fprintf(fout_case,'%s\t%s\t%s\t%s\t%s\t%s\n','N_f','Casenum',...
                'Mintime','Maxtime','numavgpoints','Rgsq_avg');
            
            
            avgcase_store = zeros(length(casearr),1); % store each average value
            
            for casecntr = 1:length(casearr) % begin case loop
                casenum = casearr(casecntr);
                
                dirname = sprintf('./../../rganalysis/outresults_dir_n_%d/%s/pdifree_%s_pdigraft_%s/Case_%d',...
                    nval,dirstr,pdifree_str,pdigraft_str,casenum);
                
                if ~exist(dirname,'dir')
                    fprintf('%s does not exist\n',dirname);
                    continue
                end
                
                fprintf('Analyzing pdi/nfree/arch/case: %g\t%d\t%s\t%d\n', pdifree,nval,dirstr,casenum);
                
                %check if file type exists
                rg_prefix = sprintf('rgavg_config_*.lammpstrj');
                rg_fylelist = dir(strcat(dirname,'/',rg_prefix));
                if min(size(rg_fylelist)) == 0
                    fprintf('No files/Empty files are found for %s\n',rg_prefix);
                    continue;
                end
                
                nfyles = numel(rg_fylelist); %number of files of the type
                reordered_ads_list = renumber_files(rg_fylelist,nfyles); % reorder file names to avoid double counting
                
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
                    
                    
                    %                     %average rg values
                    if min(data(:,1)) > set_tmin
                        
                        for minindcnt = 1:lendata %avoid double counting
                            if data(minindcnt,1) > mintstep
                                minindana = minindcnt; %minimum value at which the trajectories are separate
                                mintstep = max(data(:,1)); %new value will be the maximum value of this file
                                break;
                            end
                        end
                        sum_rgsq = sum(data(minindana:lendata,2).^2);
                        sum_across_files = sum_across_files + sum_rgsq;
                        tot_cntr_across_files = tot_cntr_across_files + length(data(minindana:lendata,2));
                        
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
                if tot_cntr_across_files ~=0
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
                rg_all(ncnt,arch_cnt)        = rg_all(ncnt,arch_cnt) + avg_for_each_casenum;
                totsamples(ncnt,arch_cnt)    = totsamples(ncnt,arch_cnt) + tot_cntr_across_files;
                avgcase_store(casecntr)      = avg_for_each_casenum;
                
                
            end % end case loop
            
            fclose(fout_case);
            
            avg_across_cases(ncnt,arch_cnt,pdi_cntr) = rg_all(ncnt,arch_cnt)/casecntr_arr(ncnt,arch_cnt);
            
            if casecntr_arr(ncnt,arch_cnt) > 1
                
                nonzerovals = avgcase_store(avgcase_store~=0);
                stderr_val = std(nonzerovals)/sqrt(length(nonzerovals));
                
            else
                
                stderr_val = -1; % Invalid for just one case
                
            end
            
            fprintf(fout_avg,'%d\t%s\t%d\t%d\t%g\t%g\n',nval,dirstr,casecntr_arr(ncnt,arch_cnt),...
                totsamples(ncnt,arch_cnt),avg_across_cases(ncnt,arch_cnt,pdi_cntr),stderr_val);
            
        end % end arch loop
        
    end % end nfree loop
    
    fclose(fout_avg);
    
    if plotrg
        h1 = figure;
        hold on
        box on
        set(gca,'FontSize',16)
        xlabel('$N_{pa}/N_{pc}$','FontSize',20,'Interpreter','Latex')
        ylabel('$\langle R_{g} \rangle^{1/2}$','FontSize',20,'Interpreter','Latex')
        
        for plcnt = 1:length(arch_arr)
            plot(nfreearr/ngraft,avg_across_cases(:,plcnt,pdi_cntr).^(0.5),'color',pclr{plcnt},...
                'Marker',msty{plcnt},'MarkerSize',8,'MarkerFaceColor',pclr{plcnt},'LineStyle','None')
            legendinfo{plcnt} = leg_arr{plcnt};
        end
        %overlay y = x line
        %         x = 0:0.5:5; y = x;
        %         plot(x,y,'LineWidth',2,'Color',orange,'LineStyle','--')
        %         ylim([min(min(avg_across_cases(:,:,pdi_cntr))) 1.2*max(max(avg_across_cases(:,:,pdi_cntr)))])
        legend(legendinfo,'FontSize',16,'Location','Best')
        legend boxoff
        saveas(h1,sprintf('./../../Figs_paper/adsorbmon_rcut_%s_pdi_%g.png',cutoff,pdifree));
        clear legendinfo
    end
    
end % end pdi free loop

fclose(fout_cons);


