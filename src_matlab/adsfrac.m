%% To compute average adsorbed fraction of chains
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
pdigraft = 1.0;
nmonfree = 30; nmongraft = 30; ngraft = 64;
cutoff = '1.50';
lz = 120; area=35^2;
set_tmax = 3e7; % maximum timestep for analysis;

%% Input flags
stddevon = 0;
adsflag  = 1;

%% Pre-calculations
rhofree = nfreearr*30/(lz*area);
pdigraft_str = num2str(pdigraft,'%1.1f');

%% Main Analysis

if adsflag % create consolidated list
    fout_cons = fopen(sprintf('./../../outfiles/overall/adsorbed_chain_consolidated_rcut_%s.dat',...
        cutoff),'w');
    fprintf(fout_cons,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n','PDI_free','N_f','Arch',...
        'Case_num','numsample_pts','min_time','max_time','avg_fraction');
end


for pdi_cntr = 1:length(pdi_freearr) % begin pdi free loop
    pdifree     = pdi_freearr(pdi_cntr);
    pdifree_str = num2str(pdifree,'%1.1f');
    
    %zero arrays for averages across cases
    casecntr_arr  = zeros(length(nfreearr),length(arch_arr));
    nadschain_all = zeros(length(nfreearr),length(arch_arr));
    totsamples    = zeros(length(nfreearr),length(arch_arr));
    
    if adsflag %create average across all cases
        fout_avg = fopen(sprintf('./../../outfiles/overall/adsorbed_chain_ave_allcases_rcut_%s_pdifree_%g.dat',...
            cutoff,pdifree),'w');
        fprintf(fout_avg,'%s\t%s\t%s\t%s\t%s\n','N_f','Arch','ncases','numsample_pts','avg_fraction');
    end
    
    for ncnt = 1:length(nfreearr) % begin nfree loop
        nval = nfreearr(ncnt);
        
        for arch_cnt = 1:length(arch_arr)  % begin arch loop
            dirstr = arch_arr{arch_cnt};
            
            dirname = sprintf('./../../sim_results/outresults_dir_n_%d/%s/pdifree_%s_pdigraft_%s',...
                nval,dirstr,pdifree_str,pdigraft_str);
            
            if ~exist(dirname,'dir')
                fprintf('%s does not exist\n',dirname);
                continue
            end
            
            if adsflag % create case-based avg outfiles
                fout_case = fopen(sprintf('./../../outfiles/n_%d/adsorbed_chain_rcut_%s_pdifree_%g_%s.dat',...
                    nval,cutoff,pdifree,dirstr),'w');
                fprintf(fout_case,'%s\t%s\t%s\t%s\t%s\t%s\n','N_f','Casenum','Mintime','Maxtime','numavgpoints','avg_fraction');
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
                    ads_prefix = sprintf('adsfracchain_rcut_%s_config_*.lammpstrj',cutoff);
                    ads_fylelist = dir(strcat(dirname,'/',ads_prefix));
                    if min(size(ads_fylelist)) == 0
                        fprintf('No files/Empty files are found for %s\n',ads_prefix);
                        continue;
                    end
                    
                    nfyles = numel(ads_fylelist); %number of files of the type
                    
                    sum_across_files = 0; tot_cntr_across_files = 0; mintime = 0; maxtime = 0;
                    for fylcnt = 1:nfyles % begin running through all files of the given type
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
                        data = importdata(ads_fylename);
                        nads_fracchain = sum(data(:,3));
                        sum_across_files = sum_across_files + nads_fracchain;
                        tot_cntr_across_files = tot_cntr_across_files + length(data(:,3));
                        
                        %find minimum and maximum time
                        if min(data(:,1)) < mintime
                            mintime = min(data(:,1));
                        end
                        if max(data(:,1)) < mintime
                            maxtime = max(data(:,1));
                        end
                        
                    end % end summing adsfrac across all files
                    
                    avg_for_each_casenum = sum_across_files/tot_cntr_across_files;
                    if maxtime > set_tmax
                        fprintf(fout_case,'%d\t%s\t%d\t%d\t%d\t%g\n',nval,casenum,...
                            mintime,maxtime,tot_cntr_across_files,avg_for_each_casenum);   
                        fprintf(fout_cons,'%g\t%d\t%s\t%d\t%d\t%d\t%d\t%g\n',pdifree,nval,...
                            dirstr,casenum,mintime,maxtime,tot_cntr_across_files,avg_for_each_casenum);
                    else
                        fprintf(fout_case,'%d\t%s\t%d\t%d\t%d\t%g\t%s\n',nval,casenum,...
                            mintime,maxtime,tot_cntr_across_files,avg_for_each_casenum,'Incomplete sampling');
                        fprintf(fout_cons,'%g\t%d\t%s\t%d\t%d\t%d\t%d\t%g\t%s\n',pdifree,nval,...
                            dirstr,casenum,mintime,maxtime,tot_cntr_across_files,avg_for_each_casenum,'Incomplete sampling');
                    end
                    
                    %save it to overall arrays
                    casecntr_arr(ncnt,arch_cnt)  = casecntr_arr(ncnt,arch_cnt) + 1;
                    nadschain_all(ncnt,arch_cnt) = nadschain_all(ncnt,arch_cnt) + avg_for_each_casenum;
                    totsamples(ncnt,arch_cnt)    = totsamples(ncnt,arch_cnt) + tot_cntr_across_files;
                    
                end %end adsorption calculation
                
            end % end case loop
            
            fclose(fout_case);
            
            avg_across_cases = nadschain_all(ncnt,arch_cnt)/casecntr_arr(ncnt,arch_cnt);
            fprintf(fout_avg,'%d\t%s\t%d\t%d\t%g\n',nval,dirstr,casecntr_arr(ncnt,arch_cnt),totsamples(ncnt,arch_cnt),avg_across_cases);
            
        end % end arch loop
        
    end % end nfree loop
    
    fclose(fout_avg);
    
end % end pdi free loop

fclose(fout_cons);