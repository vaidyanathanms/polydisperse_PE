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
nfreearr = [16;32;48;64;72;80;100];
casearr  = [1,2,3,4];
pdi_freearr = [1,1.3,1.5];
arch_arr = {'bl_bl';'bl_al';'al_bl';'al_al'};
pdigraft = 1.0;
nmonfree = 30; nmongraft = 30; ngraft = 64;
cutoff = '1.50';
lz = 120; area=35^2;

%% Input flags
stddevon = 0;
adsflag  = 1;

%% Pre-calculations
rhofree = nfreearr*30/(lz*area);
pdigraft_str = num2str(pdigraft,'%1.1f');

%% Main Analysis

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
        
        for i = length(arch_arr)  % begin arch loop
            dirstr = arch_arr{i};
            
            if adsflag % create case-based avg outfiles
                fout_case = fopen(sprintf('./../../outfiles/n_%d/adsorbed_chain_rcut_%s_pdifree_%g_%s.dat',...
                    nval,cutoff,pdifree,dirstr),'w');
                fprintf(fout_case,'%s\t%s\t%s\t%s\n','N_f','Casenum','numavgpoints','avg_fraction');
            end
            
            for casecntr = 1:length(casearr) % begin case loop
                casenum = casearr(casecntr);
                
                dirname = sprintf('./../../sim_results/ouresults_n_%d/%s/pdifree_%s_pdigraft_%s/Case_%d',...
                    nval,dirstr,pdifree_str,pdigraft_str,casenum);
                
                if ~exist(dirname,'dir')
                    fprintf('%s does not exist\n',dirname);
                    continue
                end
                
                if adsflag %begin adsorption calculation
                    
                    %check if file type exists
                    ads_prefix = sprintf('adsfracchain_rcut_%s_config_*.lammpstrj',cutoff);
                    ads_fylelist = dir(strcat(simdirname,'/',ads_prefix));
                    if min(size(ads_fylelist)) == 0
                        fprintf('No files/Empty files are found for %s\n',ads_prefix);
                        continue;
                    end
                    
                    nfyles = numel(ads_fylelist); %number of files of the type
                    
                    sum_across_files = 0; tot_cntr_across_files = 0;
                    for fylcnt = 1:nfyles % begin running through all files of the given type
                        ads_fylename = fylelist(fylcnt).name;
                        if exist(ads_fylename,'file') ~= 2
                            fprintf('%s does not exist/empty file\n',ads_fylename);
                            continue;
                        elseif struct(dir(ads_fylename)).bytes == 0
                            fprintf('Empty file: %s \n',ads_fylename);
                            continue;
                        end
                        
                        data = importdata(ads_fylename);
                        nads_fracchain = sum(data(:,3));
                        sum_across_files = sum_across_files + nads_fracchain;
                        tot_cntr_across_files = tot_cntr_across_files + length(data(:,3));
                    end % end summing adsfrac across all files
                    
                    avg_for_each_casenum = sum_across_files/tot_cntr_across_files;
                    fprintf(fout_case,'%d\t%s\t%d\t%g\n',nval,casenum,tot_cntr_across_files,avg_for_each_casenum);
                    
                    %save it to overall arrays
                    casecntr_arr(ncnt,i)  = casecntr_arr(ncnt,i) + 1;
                    nadschain_all(ncnt,i) = nadschain_all(ncnt,i) + avg_for_each_casenum;
                    totsamples(ncnt,i)    = totsamples(ncnt,i) + tot_cntr_across_files;
                    
                end %end adsorption calculation
                
            end % end case loop
    
            avg_across_cases = nadschain_all(ncnt,i)/casecntr_arr(ncnt,i);
            fprintf(fout_avg,'%d\t%s\t%d\t%d\t%g\n',nval,dirstr,casecntr_arr(ncnt,i),totsamples(ncnt,i),avg_across_cases);
            
        end % end arch loop
        
    end % end nfree loop
    
end % end pdi free loop

