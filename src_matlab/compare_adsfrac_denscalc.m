% Compare integrals of adsorbed densities and from definition that if one
% monomer is adsorbed, all the monomers are used against the calculation.

clc;
close all;
clear;
format long

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
ttestflag = 0; % write to ttest_dir the individual cases
plotads   = 0; % plot fads as a function of PDI

%% Graft details
ncharge_mons_graft = 30; % per graft chain details
ntail_mons_graft = 5;    % per graft chain details
ntot_mons_graft  = ncharge_mons_graft + ntail_mons_graft;
nch_graft = 32;

%% Zero arrays
avg_across_cases = zeros(length(nfreearr),length(arch_arr),length(pdi_freearr));

%% Pre-calculations
rhofree = nfreearr*nmonfree/(lz*area);
pdigraft_str = num2str(pdigraft,'%1.1f');
num_cases = length(casearr);
max_mw_free = 10*nmonfree; % An approximate max. Will throw error from extract_adschain.m if it is more than this value.

%% For averaging across cases (array name: avgads_molarr)
bin_wid  = 8; % THIS CAN BE DIFFERENT FROM WHAT IS IN COMPUTE_MWDIST
bin_lims = [1,max_mw_free]; % THIS IS NEEDED TO ENSURE ALL CASES ARE BINNED IN THE SAME RANGE

%% Main Analysis

s1 = create_output_dirs('./../../outfiles');
s2 = create_output_dirs('./../../outfiles/overall');

% Create consolidated list
fout_cons = fopen(sprintf('./../../monads/overall/adsorbed_mon_consolidated_rcut_%s.dat',...
    cutoff),'w');
fprintf(fout_cons,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n','PDI_free','N_f','Arch',...
    'Case_num','min_time','max_time','numsample_pts','avg_fraction');


for pdi_cntr = 1:length(pdi_freearr) % begin pdi free loop
    pdifree     = pdi_freearr(pdi_cntr);
    pdifree_str = num2str(pdifree,'%1.1f');
    
    %zero arrays for averages across cases
    casecntr_arr  = zeros(length(nfreearr),length(arch_arr));
    nadsmon_all = zeros(length(nfreearr),length(arch_arr));
    totsamples    = zeros(length(nfreearr),length(arch_arr));
    
    % Create average across all cases
    fout_avg = fopen(sprintf('./../../monads/overall/adsorbedmon_wrtch_ave_allcases_rcut_%s_pdifree_%g.dat',...
        cutoff,pdifree),'w');
    fprintf(fout_avg,'%s\t%s\t%s\t%s\t%s\t%s\n','N_f','Arch','ncases','numsample_pts',...
        'avg_fraction','StdErrMean');
    
    % Create case-based avg outfiles for SI data
    fout_sidata = fopen(sprintf('./../../monads/overall/sidata_adsorbedmon_wrtch_rcut_%s_pdifree_%g.dat',...
        cutoff,pdifree),'w');
    fprintf(fout_sidata,'%s\t%s\t%s\t%s\t%s\n','\DJ$_{\rm{ideal}}$','$N_{pa}$','Arch',...
        'Case \#','avg_fraction');
    
    
    for ncnt = 1:length(nfreearr) % begin nfree loop
        nval = nfreearr(ncnt);
        s2 = create_output_dirs(sprintf('./../../outfiles/overall/n_%d',nval));

        for arch_cnt = 1:length(arch_arr)  % begin arch loop
            dirstr = arch_arr{arch_cnt};
            
            dirname = sprintf('./../../sim_results/outresults_dir_n_%d/%s/pdifree_%s_pdigraft_%s',...
                nval,dirstr,pdifree_str,pdigraft_str);
            
            if ~exist(dirname,'dir')
                fprintf('%s does not exist\n',dirname);
                continue
            end
            
            s3 = create_output_dirs(sprintf('./../../monads/n_%d',nval));
            % Create case-based avg outfiles
            fout_case = fopen(sprintf('./../../monads/n_%d/adsorbedmon_wrtch_rcut_%s_pdifree_%g_%s.dat',...
                nval,cutoff,pdifree,dirstr),'w');
            fprintf(fout_case,'%s\t%s\t%s\t%s\t%s\t%s\n','N_f','Casenum',...
                'Mintime','Maxtime','numavgpoints','avg_fraction');
            
            
            if ttestflag %write all cases into separate folders in ttest_dir
                
                fout_ttest = fopen(sprintf('./../../monads/ttest_dir/n_%d/adsfracmon_chmw_rcut_%s_pdifree_%g_arch_%s.dat',...
                    nval,cutoff,pdifree,dirstr),'w');
                fprintf(fout_ttest,'%s\n','avg_fraction: #of columns correspond to the number of cases');
                
            end
            
            
            avgcase_store = zeros(length(casearr),1); % store each average value
            
            for casecntr = 1:length(casearr) % begin case loop
                casenum = casearr(casecntr);
                
                dirname = sprintf('./../../sim_results/outresults_dir_n_%d/%s/pdifree_%s_pdigraft_%s/Case_%d',...
                    nval,dirstr,pdifree_str,pdigraft_str,casenum);
                
                if ~exist(dirname,'dir')
                    fprintf('%s does not exist\n',dirname);
                    continue
                end
                
                fprintf('Analyzing pdi/nfree/arch/case: %g\t%d\t%s\t%d\n', pdifree,nval,dirstr,casenum);
                
                %check if file type exists
                ads_prefix = sprintf('adsfrac_chmw_rcut_%s_config_*.lammpstrj',cutoff);
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
                        nads_fracmon = sum(data(minindana:lendata,2));
                        sum_across_files = sum_across_files + nads_fracmon;
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
                nadsmon_all(ncnt,arch_cnt)   = nadsmon_all(ncnt,arch_cnt) + avg_for_each_casenum;
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
                avg_across_cases(ncnt,arch_cnt,pdi_cntr) = nadsmon_all(ncnt,arch_cnt)/length(nonzerovals);

            else
                
                stderr_val = -1; % Invalid for just one case
                
            end
            
            fprintf(fout_avg,'%d\t%s\t%d\t%d\t%g\t%g\n',nval,dirstr,casecntr_arr(ncnt,arch_cnt),totsamples(ncnt,arch_cnt),avg_across_cases(ncnt,arch_cnt,pdi_cntr),stderr_val);
            
        end % end arch loop
        
            end % end nfree loop
    
    fclose(fout_avg);
    
end
                    dirname = sprintf('./../../sim_results/outresults_dir_n_%d/%s/pdifree_%s_pdigraft_%s/Case_%d',...
                        nplot,dirstr,pdifree_str,pdigraft_str,casenum);
                    if ~exist(dirname,'dir')
                        fprintf('%s does not exist\n',dirname);
                        continue
                    end
                    rho_prefix = 'densadsch_config_*.lammpstrj';
                    rho_fylelist = dir(strcat(dirname,'/',rho_prefix));
                    if min(size(rho_fylelist)) == 0
                        fprintf('No files/Empty files are found for %s\n',rho_prefix);
                        continue;
                    end
                    
                    nfyles = numel(rho_fylelist); %number of files of the type
                    if nfyles == 0
                        fprintf('Did not find files of the type grpdens_config* in %d\t%s\n', nplot,dirstr);
                        continue;
                    else
                        [latest_fyleindex] = find_latest_fyle(rho_fylelist,nfyles);
                    end
                    
                    rho_fylename = strcat(dirname,'/',rho_fylelist(latest_fyleindex).name);
                    if exist(rho_fylename,'file') ~= 2
                        fprintf('%s does not exist/empty file\n',rho_fylename);
                        continue;
                    elseif struct(dir(rho_fylename)).bytes == 0
                        fprintf('Empty file: %s \n',rho_fylename);
                        continue;
                    end
                    fprintf('Plotting density profiles using %s for n_pa = %d\n', rho_fylename, nplot);
                    
                    all_data = importdata(rho_fylename,' ',1);
                    fld = all_data.data;
                    rdata = fld(:,1); lz = fld(1,1) + fld(length(fld(:,1)),1); %because of binning
                    nbins = length(rdata(:,1));
                    
                    % sanity check for length of data
                    if casecntr == 1
                        nbins_old = nbins;
                    else
                        if nbins_old ~= nbins
                            fprintf('ERROR: The values of bins do not match %s\t%d\t%g\n',dirstr,nplot,pdifree)
                            error('CHECK nbins')
                        else
                            nbins_old = nbins;
                        end
                    end
                    
                    avg_rho_neutral = fld(:,2) + avg_rho_neutral;
                    avg_rho_charged  = fld(:,3) + avg_rho_charged;
                    
                    totcases = totcases + 1;
                    
                end
                
                avg_rho_neutral =  avg_rho_neutral/totcases;
                avg_rho_charged  =  avg_rho_charged/totcases;