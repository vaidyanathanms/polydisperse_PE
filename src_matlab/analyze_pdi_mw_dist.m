%% To compute pdi, initial and final molecular weight distributions
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
nfreearr = [16;32;64;96;128;150];
casearr  = [1,2,3,4];
pdi_freearr = [1,1.3,1.5];
arch_arr = {'bl_bl','bl_al','al_bl','al_al'};
leg_arr  = {'Block-Block','Block-Alter','Alter-Block','Alter-Alter'}; % ALWAYS CHECK for correspondence with arch_arr
pdigraft = 1.0;
nmonfree = 30; nmongraft = 30; ngraft = 32;
cutoff = '1.50';
lz = 120; area=35^2;
set_tmax = 3e7; % maximum timestep for analysis;
nfreemons = 30;

%% Input flags
pdiflag  = 1;
mwdflag  = 1;

%% Zero arrays
avg_across_cases = zeros(length(nfreearr),length(arch_arr),length(pdi_freearr));

%% Pre-calculations
rhofree = nfreearr*nfreemons/(lz*area);
pdigraft_str = num2str(pdigraft,'%1.1f');

%% Main Analysis

if pdiflag % create consolidated list
    fout_cons = fopen(sprintf('./../../outfiles/overall/pdi_%s.dat',...
        cutoff),'w');
    fprintf(fout_cons,'%s\t%s\t%s\t%s\t%s\t%s\n','N_f','Arch','Case_num',...
        'Mw_free','Mn_free','PDI_free');
end

for pdi_cntr = 1:length(pdi_freearr) % begin pdi free loop
    pdifree     = pdi_freearr(pdi_cntr);
    pdifree_str = num2str(pdifree,'%1.1f');
    
    %zero arrays for averages across cases
    casecntr_arr  = zeros(length(nfreearr),length(arch_arr));
    nadschain_all = zeros(length(nfreearr),length(arch_arr));
    totsamples    = zeros(length(nfreearr),length(arch_arr));
    
    if pdiflag %create average across all cases
        fout_avg = fopen(sprintf('./../../outfiles/overall/pdi_ave_allcases_rcut_%s.dat',...
            cutoff),'w');
        fprintf(fout_avg,'%s\t%s\t%s\t%s\n','N_f','Arch','ncases','pdi_avg');
    end
    
    for ncnt = 1:length(nfreearr) % begin nfree loop
        nval = nfreearr(ncnt);
        
        for arch_cnt = 1:length(arch_arr)  % begin arch loop
            dirstr = arch_arr{arch_cnt};
            
            dirname = sprintf('./../../data_all_dir/n_%d/%s/pdifree_%s_pdigraft_%s',...
                nval,dirstr,pdifree_str,pdigraft_str);
            
            if ~exist(dirname,'dir')
                fprintf('%s does not exist\n',dirname);
                continue
            end
            
            for casecntr = 1:length(casearr) % begin case loop
                casenum = casearr(casecntr);
                
                dirname = sprintf('./../../sim_results/outresults_dir_n_%d/%s/pdifree_%s_pdigraft_%s/Case_%d',...
                    nval,dirstr,pdifree_str,pdigraft_str,casenum);
                
                if ~exist(dirname,'dir')
                    fprintf('%s does not exist\n',dirname);
                    continue
                end
                
                fprintf('Analyzing pdi/nfree/arch/case: %g\t%d\t%s\t%d\n', pdifree,nval,dirstr,casenum);
                if adsflag % begin adsorption calculation
                    
                    %check if file type exists
                    ads_prefix = sprintf('PEinitdata');
                    ads_fylelist = dir(strcat(dirname,'/',ads_prefix));
                    if min(size(ads_fylelist)) == 0
                        fprintf('No files/Empty files are found for %s\n',ads_prefix);
                        continue;
                    end
                    
                    nfyles = numel(ads_fylelist); %number of files of the type
                    
                    if nfyles ~= 1
                        fprintf('WARNING: Found %d initial datafiles', nfyles);
                        continue;
                    end
                    
                    sum_across_files = 0; tot_cntr_across_files = 0; mintime = 10^10; maxtime = 0;
                    % begin running through the datafile
                    ads_fylename = strcat(dirname,'/',ads_fylelist(fylcnt).name);
                    if exist(ads_fylename,'file') ~= 2
                        fprintf('%s does not exist/empty file\n',ads_fylename);
                        continue;
                    elseif struct(dir(ads_fylename)).bytes == 0
                        fprintf('Empty file: %s \n',ads_fylename);
                        continue;
                    end
                    
                    %average adsorption values
                    fprintf('Analyzing %s\n', ads_fylename);
                    molarr_cnt = analyze_datafile(ads_fylename,nval);
                    %compute_init_pdi();
                    %compute_mwdist();
                    
                    %save it to overall arrays
                    %casecntr_arr(ncnt,arch_cnt)  = casecntr_arr(ncnt,arch_cnt) + 1;
                    %nadschain_all(ncnt,arch_cnt) = nadschain_all(ncnt,arch_cnt) + avg_for_each_casenum;
                    %totsamples(ncnt,arch_cnt)    = totsamples(ncnt,arch_cnt) + tot_cntr_across_files;
                    
                end %end adsorption calculation
                
            end % end case loop
            
            fclose(fout_case);
            
            avg_across_cases(ncnt,arch_cnt,pdi_cntr) = nadschain_all(ncnt,arch_cnt)/casecntr_arr(ncnt,arch_cnt);
            fprintf(fout_avg,'%d\t%s\t%d\t%d\t%g\n',nval,dirstr,casecntr_arr(ncnt,arch_cnt),totsamples(ncnt,arch_cnt),avg_across_cases(ncnt,arch_cnt,pdi_cntr));
            
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
        saveas(h1,sprintf('./../../all_figures/adsorbchain_rcut_%s_pdi_%g.png',cutoff,pdifree));
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
        saveas(h1,sprintf('./../../all_figures/adsorbchain_rcut_%s_nval_%g.png',cutoff,nval));
        clear legendinfo
        
    end
    
end
