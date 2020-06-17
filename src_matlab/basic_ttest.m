%% To compare between architectures and within the same architecture for correlation

% Added two sample unequal variance t-test
% Ref: https://www.mathworks.com/help/stats/ttest2.html#btrj_js-1

% Need to add one sample equal variance t-test

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
cutoff_arr = {'1.30','1.50'};
lz = 120; area=35^2;
set_tmax = 3e7; % maximum timestep for analysis;

%% Input flags
ttest_1_flag = 0;
ttest_2_flag = 0;

%% Zero arrays
avg_across_cases = zeros(length(nfreearr),length(arch_arr),length(pdi_freearr));
num_of_cases     = zeros(length(nfreearr),length(arch_arr),length(pdi_freearr));
nval_from_fyle   = zeros(length(nfreearr),length(arch_arr),length(pdi_freearr));
pdi_from_fyle    = zeros(length(nfreearr),length(arch_arr),length(pdi_freearr));

%% Pre-calculations
rhofree = nfreearr*30/(lz*area);
pdigraft_str = num2str(pdigraft,'%1.1f');
err_tol = 1e-10; % for finding elements in a real array

%% Main Analysis

for rcutcntr = 1:length(rcut_arr) % begin rcut loop
    cutoff = cutoff_arr{rcutcntr};
    
    % read main file
    fylename = sprintf(sprintf('./../../outfiles/overall/adsorbed_chain_consolidated_rcut_%s.dat',...
        cutoff));
    fin_cons = fopen(fylename,'r');
    
    if fin_cons <= 0 % check for average list
        fprintf('%s does not exist', fylename);
        continue;
    end
    
    % read main file
    fylename = sprintf(sprintf('./../../outfiles/overall/cross_check_adsorbed_chain_consolidated_rcut_%s.dat',...
        cutoff));
    fout_vals = fopen(fylename,'w');
    
    
    headerflag = -1; % find headerflag
    pdi_col = -1; nf_col = -1; arch_col = -1; casenum_col = -1;
    numsam_col = -1; avgf_col = -1;
    
    while true
        
        if headerflag == -1
            
            tline = fgetl(fout_cons);
            if ~ischar(tline) || isempty(tline)
                continue;
            end
            
            spl_tline = strsplit(strtrim(tline));
            len_tline = length(spl_tline);
            
            if len_tline ~= 8
                fprintf('Unknown number of keywords in %s \n', tline);
                continue;
            end
            
            for word_cnt = 1:len_tline % add more if keywords are added
                
                if strcmp(spl_tline{word_cnt},'PDI_free') % for pdi_free keyword
                    pdi_col = word_cnt;
                elseif strcmp(spl_tline{word_cnt},'N_f') % for N_f keyword
                    nf_col = word_cnt;
                elseif strcmp(spl_tline{word_cnt},'Arch') % for Arch keyword
                    arch_col = word_cnt;
                elseif strcmp(spl_tline{word_cnt},'Case_num') % for Case_num keyword
                    casenum_col = word_cnt;
                elseif strcmp(spl_tline{word_cnt},'numsample_pts') % for number of samples keyword
                    numsam_col = word_cnt;
                elseif strcmp(spl_tline{word_cnt},'avg_fraction') % for avg_fraction keyword
                    avgf_col = word_cnt;
                end
                
            end
            
            headerflag = 1; % found headerflag
            
        else
            
            tline = fgetl(fout_cons);
            spl_tline = strsplit(strtrim(tline));
            if ~ischar(tline) || isempty(tline)
                continue;
            end
            
            spl_tline = strsplit(strtrim(tline));
            len_tline = length(spl_tline);
            
            if len_tline ~= 8 && len_tline ~= 9
                fprintf('Unknown number of keywords in %s \n', tline);
                continue;
            end
            
            % parse elements
            nval   = str2double(spl_tline{nf_col});
            dirstr = spl_tline{arch_col};
            pdival = str2double(spl_tline{pdi_col});
            
            % cross check with the main array
            if ismember(nval, nfreearr)
                nval_index = find(nval==nfreearr);
            else
                fprintf('did not find %d in nfreearr\n', nval_index)
                continue;
            end
            
            if contains(arch_arr,dirstr)
                arch_index = find(contains(arch_arr,dirstr));
            else
                fprintf('did not find %s in arch_arr\n', dirstr)
                continue;
            end
            
            if ismembertol(pdival,pdi_freearr,err_tol)
                pdi_index = find(pdival,pdi_freearr);
            else
               fprintf('did not find %g in pdifree within a tolerance of %g\n', pdival,err_tol)
               continue;
            end 
            
            avg_across_cases(nval_index,arch_index,pdi_index) = str2double(spl_tline{avgf_col});
            nval_from_fyle(nval_index,arch_index,pdi_index) = nval;
            pdi_from_fyle(nval_index,arch_index,pdi_index) = pdival;
            num_of_cases(nval_index,arch_index,pdi_index) = num_cases(nval_index,arch_index,pdi_index) + 1;
            
            fin_vals('%g\t%d\t%s\t%s\t%g\n',pdi_from_fyle(nval_index,arch_index,pdi_index),...
                nval_from_fyle(nval_index,arch_index,pdi_index),dirstr,...
                num_of_cases(nval_index,arch_index,pdi_index)
                
            
        end % end if loop for reading lines
            
    end % end file-read loop
    
    fclose(fin_cons);
    
end % end rcut loop
