%% To plot output MW distribution
% V2.0 changed binning algorithm, see NOTE 2 and NOTE 3
% Results after NOTE 2 will mostly be irrelevant (only for cross-checking)

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
nch_freearr = [32]%;64;128;150];
casearr  = [1]%;2;3;4];
pdi_freearr = [1.5];
arch_arr  = {'bl_bl'}%;'bl_al';'al_bl';'al_al'};
leg_arr   = {'Block-Block'}%;'Block-Alter';'Alter-Block';'Alter-Alter'}; % ALWAYS CHECK for correspondence with arch_arr for legends
pdigraft  = 1.0;
nfreemons = 30;
ngraft_ch = 32; % Number of graft chains
cutoff = '1.50';
lz = 120;
area = 35^2;
set_tmax = 3e7; % maximum timestep for analysis;

%% Graft details

ncharge_mons_graft = 30; % per graft chain details
ntail_mons_graft = 5;    % per graft chain details
ntot_mons_graft  = ncharge_mons_graft + ntail_mons_graft;
nch_graft = 32;

%% Zero arrays
avg_across_cases = zeros(length(nch_freearr),length(arch_arr),length(pdi_freearr));

%% Pre-calculations
rhofree = nch_freearr*nfreemons/(lz*area);
pdigraft_str = num2str(pdigraft,'%1.1f');
num_cases = length(casearr);
max_mw_free = 10*nfreemons; % An approximate max. Will throw error from extract_adschain.m if it is more than this value.

%% For averaging across cases
bin_wid  = 8; % THIS CAN BE DIFFERENT FROM WHAT IS IN COMPUTE_MWDIST
bin_lims = [1,max_mw_free]; % THIS IS NEEDED TO ENSURE ALL CASES ARE BINNED IN THE SAME RANGE
bin_edges = []

%% Main Analysis

for ncnt = 1:length(nch_freearr) % begin nfree loop
    nval = nch_freearr(ncnt);
    
    for pdi_cntr = 1:length(pdi_freearr) % begin pdi free loop
        ref_pdifree     = pdi_freearr(pdi_cntr);
        pdifree_str     = num2str(ref_pdifree,'%1.1f');
        
        for arch_cnt = 1:length(arch_arr)  % begin arch loop
            dirstr = arch_arr{arch_cnt};
            
            dirname = sprintf('./../../data_all_dir/n_%d/%s/pdifree%s_pdigraft_%s',...
                nval,dirstr,pdifree_str,pdigraft_str);
            
            if ~exist(dirname,'dir')
                fprintf('%s does not exist\n',dirname);
                continue
            end
            
            % Output MWD averaged across cases
            % avgads_molarr consists of the total number a chain of a given
            % MW is adsorbed and is totalled across all the files for a
            % given case and across cases for a given graft-free
            % architecure.
            avgads_molarr = zeros(max_mw_free,2); % to compute average distribution; maximum size of the array should be equal to the expected max MW
            avgads_molarr(:,1) = 1:max_mw_free; init_index_avgads = 1;
            avg_of_all_ads_mw_arr = zeros(1000,1); % The number 1000 is by default. Will weed zeros at the end. This will append all the MWs of the adsorbed chains
            favg_dist = fopen(sprintf('./../../outfiles/overall/out_mwdist_n_%d_pdi_%g_%s_rcut_%s.dat',...
                nval,ref_pdifree,dirstr,cutoff),'w');
            nframes_arch = 0; % Total frames per arch: sum across different cases and different files.
            
            % Avg input MW for normalization: unlike avgads the size here
            % is fixed
            all_INIT_mw_arr = zeros(length(casearr)*nval,1); %
            
            for casecntr = 1:length(casearr) % begin case loop
                casenum = casearr(casecntr);
                
                fprintf('Analyzing pdi/nfree/arch/case: %g\t%d\t%s\t%d\n', ref_pdifree,nval,dirstr,casenum);
                
                %read molecular details from input datafile and remap molecular IDs of free chains
                inp_fylename = sprintf('./../../data_all_dir/n_%d/%s/pdifree%s_pdigraft_%s/Case_%d/PEinitdata.txt',...
                    nval,dirstr,pdifree_str,pdigraft_str,casenum);
                if exist(inp_fylename,'file') ~= 2
                    fprintf('%s does not exist/empty file\n',inp_fylename);
                    continue;
                end
                % molarr outputs the remapped ID of chain, actual ID and MW
                % of each chain
                molarr = analyze_datafile(inp_fylename,nval,nch_graft); % extract molecular details for comparison at the end
                start_index = (casenum-1)*nval + 1; fin_index = casenum*nval;
                all_INIT_mw_arr(start_index:fin_index) = molarr(:,3); % contains the initial MW of all free chains
                
                % Now start analyzing all adsorbed chain files
                dirname = sprintf('./../../sim_results/outresults_dir_n_%d/%s/pdifree_%s_pdigraft_%s/Case_%d',...
                    nval,dirstr,pdifree_str,pdigraft_str,casenum);
                
                if ~exist(dirname,'dir')
                    fprintf('%s does not exist\n',dirname);
                    continue
                end
                
                % check if adsorbed chain file type exists
                ads_prefix = sprintf('chainadsval_rcut_%s_config_*.lammpstrj',cutoff);
                ads_fylelist = dir(strcat(dirname,'/',ads_prefix));
                if min(size(ads_fylelist)) == 0
                    fprintf('No files/Empty files are found for %s\n',ads_prefix);
                    continue;
                end
                
                nfyles = numel(ads_fylelist); %number of files of the type
                nframes_case = 0; % total frames per case: sum across all files.
                
                % begin running through the chainadsfile
                for fylcnt = 1:nfyles
                    ads_fylename = strcat(dirname,'/',ads_fylelist(fylcnt).name);
                    if exist(ads_fylename,'file') ~= 2
                        fprintf('%s does not exist/empty file\n',ads_fylename);
                        continue;
                    elseif struct(dir(ads_fylename)).bytes == 0
                        fprintf('Empty file: %s \n',ads_fylename);
                        continue;
                    end
                    
                    % extract adsorbed chain details: o/p contains
                    % frequency of all adsorbed chains, MW of each adsorbed
                    % chain and number of frames analyzed
                    [adsfreechains_arr,all_ads_mw_arr,num_frames] = extract_adschain(ads_fylename,max_mw_free);
                    nframes_case = nframes_case + num_frames;
                    
                    % copy all data into an average array
                    len_ads_arr = length(all_ads_mw_arr);
                    fin_index_avgads = init_index_avg_ads + len_ads_arr - 1;
                    avg_of_all_ads_mw_arr(init_index_avgads:fin_index_avgads) = all_ads_mw_arr(:,1); % contains the MW of all adsorbed chains
                    init_index_avg_ads = init_index_avg_ads + len_ads_arr;
                    
                    %For compute_mwdist(), the size of the array will be
                    %equal to the maximum MW and the IDs corresponding to
                    %the chains do not matter -- especially while taking an
                    %average, since for each case, the distribution of MWs
                    %will be different.
                    
                    %NOTE 2: For individual cases, Output is not added with respect to anything.
                    %If the output distribution is compared to the input,
                    %the height of the curve will show the difference
                    %between the input and output distributions.
                    outdist = compute_mwdist(adsfreechains_arr,2); % compute and write individual distribution
                    distoutfyle = strcat(dirname,'/','ads_distout_details.dat');
                    fdist = fopen(distoutfyle,'w');
                    left_edges = outdist.BinEdges(1:outdist.NumBins);
                    right_edges = outdist.BinEdges(2:outdist.NumBins+1);
                    all_values  = outdist.Values;
                    fprintf(fdist,'%s\t%d\n','NumBins',outdist.NumBins);
                    fprintf(fdist,'%s\t%s\t%s\n','leftEdge','RightEdge', 'Normalized Value');
                    fprintf(fdist,'%d\t%d\t%d\n',[left_edges' right_edges' outdist.Values']');
                    fclose(fdist);
                    
                    
                end % end mol. wt distribution calculation
                
                nframes_arch = nframes_arch + nframes_case;
                
            end % end case loop
            
            % Store into avg arrays
            % NOTE 3: Easiest way is to add all the adsorbed mol wts and
            % then count them after sorting. Subsequently divide by the
            % number of times they occur in the initial configuration
            
            [avgoutdist,avgprob] = find_distribution_of_mw(avg_of_all_ads_mw_arr(:,1),all_INIT_mw_arr(:,1),nframes_arch,bin_wid);
            
            
            % find avg probability of adsorption
            fclose(favg_dist);
            
        end % end arch loop
        
    end % end nfree loop
    
end % end pdi free loop