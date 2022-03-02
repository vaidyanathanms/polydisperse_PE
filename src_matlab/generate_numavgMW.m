%% To GENERATE number averaged MW data as a function of time from the adsorbed chain data and the adsorbed MW data
% Definition: If any of the monomers in a chain is near the prescribed
% distance, the number of monomers adsorbed = the length of the chain.
% Run compute_numavgMW.m for plotting the averages after this
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
nfreearr = [128];
casearr  = [4];
pdi_freearr = [1];
arch_arr = {'bl_bl','al_al'};
leg_arr  = {'Block-Block','Alter-Alter'}; % ALWAYS CHECK for correspondence with arch_arr
pdigraft = 1.0;
nmonfree = 30; nmongraft = 30; ngraft = 32;
cutoff = '1.50';
lz = 120; area=35^2;

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

%% Main Analysis

s1 = create_output_dirs('./../../outfiles');
s2 = create_output_dirs('./../../outfiles/overall');

for pdi_cntr = 1:length(pdi_freearr) % begin pdi free loop
    pdifree     = pdi_freearr(pdi_cntr);
    pdifree_str = num2str(pdifree,'%1.1f');    
    
    for ncnt = 1:length(nfreearr) % begin nfree loop
        nval = nfreearr(ncnt);
        s2 = create_output_dirs(sprintf('./../../outfiles/overall/n_%d',nval));

        for arch_cnt = 1:length(arch_arr)  % begin arch loop
            dirstr = arch_arr{arch_cnt};
            
            dirname = sprintf('./../../numavg_mw/outresults_dir_n_%d/%s/pdifree_%s_pdigraft_%s',...
                nval,dirstr,pdifree_str,pdigraft_str);
            
            if ~exist(dirname,'dir')
                fprintf('%s does not exist\n',dirname);
                continue
            end
            
            for casecntr = 1:length(casearr) % begin case loop
                casenum = casearr(casecntr);
                
                dirname = sprintf('./../../numavg_mw/outresults_dir_n_%d/%s/pdifree_%s_pdigraft_%s/Case_%d',...
                    nval,dirstr,pdifree_str,pdigraft_str,casenum);
                
                if ~exist(dirname,'dir')
                    fprintf('%s does not exist\n',dirname);
                    continue
                end
                
                fprintf('Analyzing pdi/nfree/arch/case: %g\t%d\t%s\t%d\n', pdifree,nval,dirstr,casenum);
                
                %check if file with MW adsorbed data exists
                ads_prefix = sprintf('adsfrac_chmw_rcut_%s_config_*.lammpstrj',cutoff);
                ads_fylelist_MW = dir(strcat(dirname,'/',ads_prefix));
                if min(size(ads_fylelist_MW)) == 0
                    fprintf('No files/Empty files are found for %s\n',ads_prefix);
                    continue;
                end
                
                nfyles_MW = numel(ads_fylelist_MW); %number of files of the type
                reordered_ads_list_MW = renumber_files(ads_fylelist_MW,nfyles_MW); % reorder file names to avoid double counting
                
                %check if file with average fraction of chains adsorbed exist
                ads_prefix = sprintf('adsfracchain_rcut_%s_config_*.lammpstrj',cutoff);
                ads_fylelist_ch = dir(strcat(dirname,'/',ads_prefix));
                if min(size(ads_fylelist_ch)) == 0
                    fprintf('No files/Empty files are found for %s\n',ads_prefix);
                    continue;
                end
                
                nfyles_ch = numel(ads_fylelist_ch); %number of files of the type
                reordered_ads_list_ch = renumber_files(ads_fylelist_ch,nfyles_ch); % reorder file names to avoid double counting
                
                if nfyles_MW ~= nfyles_ch
                    fprintf('Mismatch in number of files between adsorbed MW and adsorbed chains %d\t%d\n',nfyles_MW,nfyles_ch);
                    error('Check files');
                else
                    nfyles = nfyles_MW;
                end
                
                sum_across_files = 0; tot_cntr_across_files = 0; mintime = 10^10; maxtime = 0;
                mintstep = 0; 
                
                for fylcnt = 1:nfyles % begin running through all files of the given type
                    timeval_MW = extract_timeval(reordered_ads_list_MW{fylcnt});
                    timeval_ch = extract_timeval(reordered_ads_list_MW{fylcnt});
                    
                    if timeval_MW ~= timeval_ch
                        fprintf('Mismatch in timevalue%s\t%s\n',timeval_MW,timeval_ch);
                        error('Check files');
                    else
                        timeval = timeval_MW;
                    end
                    
                    adsout = sprintf('adsfrac_numavg_rcut_%s_config_%d.lammpstrj',cutoff,timeval);
                    ftimedata = fopen(strcat(dirname,'/',adsout),'w');
                    fprintf(ftimedata,'%s\t%s\t%s\t%s\n','Time','fadsMW','fadsch','fadsNumAvgMW');
                    
                    
                    ads_fylename_MW = strcat(dirname,'/',reordered_ads_list_MW{fylcnt});
                    
                    if exist(ads_fylename_MW,'file') ~= 2
                        fprintf('%s does not exist/empty file\n',ads_fylename_MW);
                        continue;
                    elseif struct(dir(ads_fylename_MW)).bytes == 0
                        fprintf('Empty file: %s \n',ads_fylename_MW);
                        continue;
                    end
                   
                    ads_fylename_ch = strcat(dirname,'/',reordered_ads_list_ch{fylcnt});
                    
                    if exist(ads_fylename_ch,'file') ~= 2
                        fprintf('%s does not exist/empty file\n',ads_fylename_ch);
                        continue;
                    elseif struct(dir(ads_fylename_ch)).bytes == 0
                        fprintf('Empty file: %s \n',ads_fylename_ch);
                        continue;
                    end
                    
                    fprintf('Analyzing %s\t%s\n', ads_fylename_MW,ads_fylename_ch);
                    data_MW = importdata(ads_fylename_MW);
                    data_ch = importdata(ads_fylename_ch);
                    
                    lendata_MW = length(data_MW(:,1)); lendata_ch = length(data_ch(:,1));
                    if lendata_MW ~= lendata_ch
                        lendata = min(lendata_MW,lendata_ch); % salvaging what I have
                    else
                        lendata = lendata_MW;
                    end

                    for ratcnt = 1:lendata
                        if data_ch(ratcnt,1) ~= data_MW(ratcnt,1)
                            fprintf('Mismatch in timestep %d\t%d\n',data_ch(ratcnt,1), data_MW(ratcnt,1))
                            error('Check files');
                        end
                        nads_ch  = data_ch(ratcnt,2);
                        nads_MW  = data_MW(ratcnt,2);
                        nads_rat = nads_MW/nads_ch;
                        fprintf(ftimedata,'%d\t%d\t%d\t%g\n',data_ch(ratcnt,1),nads_MW,nads_ch,nads_rat);
                    end
                         
                end % end analyzing for each file in case directory
                
                fclose(ftimedata);
                
            end % end case loop
                        
        end % end arch loop
        
    end % end nfree loop
        
end % end pdi free loop
