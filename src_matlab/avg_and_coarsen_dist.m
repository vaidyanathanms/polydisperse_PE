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

%% Flags
avg_flag = 1;
plt_flag = 1;

%% Input data
nch_freearr = [32]%;64;128;150];
casearr  = [1;2;3;4];
pdi_freearr = [1.5];
arch_arr  = {'bl_bl'}%;'bl_al';'al_bl';'al_al'};
leg_arr   = {'Block-Block'}%;'Block-Alter';'Alter-Block';'Alter-Alter'}; % ALWAYS CHECK for correspondence with arch_arr for legends
pdigraft  = 1.0;
nfreemons = 30;
ngraft_ch = 32; % Number of graft chains
cutoff = '1.50';
lz = 120;
area = 35^2;

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
                
                finp_data = fopen(sprintf('./../../outfiles/overall/out_mwdist_n_%d_pdi_%g_%s_rcut_%s.dat',...
                    nval,ref_pdifree,dirstr,cutoff),'r');
                if finp_data <= 0
                    fprintf('ERROR: %s not found\n', finp_data);
                    continue;
                end
                
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
                    tline = fgetl(finp_data);
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
                        else % if it is already present, add extra flag to third column and the normalized value to second column so that it can be divided at the end
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


%% Compute coarsened distribution



