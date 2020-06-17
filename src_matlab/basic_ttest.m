%% To compare between architectures and within the same architecture for correlation

% V2.0: Use ref1 and ref2 (fix n_f and pdi_free and check for different
% correlations). Read from ttest_dir
% Added one/two sample unequal variance t-test
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
nfreearr = [16]%;32;64;96;128;150];
casearr  = [1,2,3,4];
pdi_freearr = [1.5];
ref_arch_arr1 = {'bl_bl','bl_al','al_bl','al_al'};
ref_arch_arr2 = {'bl_bl','bl_al','al_bl','al_al'};
pdigraft = 1.0;
nmonfree = 30; nmongraft = 30; ngraft = 32;
cutoff_arr = {'1.30','1.50'};
lz = 120; area=35^2;
set_tmax = 3e7; % maximum timestep for analysis;

%% Zero arrays
avg_across_cases = zeros(length(nfreearr),length(ref_arch_arr1),length(pdi_freearr));
num_of_cases     = zeros(length(nfreearr),length(ref_arch_arr1),length(pdi_freearr));
nval_from_fyle   = zeros(length(nfreearr),length(ref_arch_arr1),length(pdi_freearr));
pdi_from_fyle    = zeros(length(nfreearr),length(ref_arch_arr1),length(pdi_freearr));

%% Pre-calculations
rhofree = nfreearr*30/(lz*area);
pdigraft_str = num2str(pdigraft,'%1.1f');
err_tol = 1e-10; % for finding elements in a real array

%% Main Analysis

for ncntr = 1:length(nfree_arr)
    nval = nfree_arr(ncntr);
    % write t-test file
    fylename = sprintf(sprintf('./../../ttest_dir/n_%d/ttestvals_rcut_%s.dat',...
        cutoff));
    fout_vals = fopen(fylename,'w');
    
    %out format:
    %pdival arch1 c1 .. c4 cavg tab tab arch2 tab tab arch2 c1 .. c4 cavg
    %tab tab tvalue
    
    fprintf('%s\t%s\t%s\t%s\t%s\t%s\t%s\t\t%s\t%s\t%s\t%s\t%s\t%s\t\t%s\n',...
        'pdival','arch1','c1','c2','c3','c4','cavg',...
        'arch2','c1','c2','c3','c4','cavg','tvalue')
    
    for pdicntr = 1:length(pdi_freearr)
        pdival = pdi_freearr(pdicntr);
        
        for rcutcntr = 1:length(rcut_arr) % begin rcut loop
            cutoff = cutoff_arr{rcutcntr};

            for arch_cntr_1 = 1:length(ref_arch_arr1) % begin arch loop

                for arch_cntr_2 = 1:length(ref_arch_arr1)
                    % read file written by adsfrac
                    fylename = sprintf(sprintf('./../../ttest_dir/n_%d/adsfrac_rcut_%s_pdifree_%g_arch_%s.dat',...
                        nval,cutoff,pdifree,dirstr));
                    fin_main = fopen(fylename,'r');
                    
                    if fin_main <= 0 % check for average list
                        fprintf('%s does not exist', fylename);
                        continue;
                    end
                    
                    fgetl(fin_main); %skip first line
                    tline = fgetl(fin_main);
                    spl_tline = strsplit(strtrim(tline));
                    len_tline = length(spl_tline);
                    
                    avg_fvals = zeros(len_tline,1); % corresponding to num of cases in the file.
                    
                    for colcntr = 1:len_tline
                        avg_fvals(col_cntr)   = str2double(spl_tline{colcntr});
                    end
                    
                end    
                
                
            
            
            
            
        end % end file-read loop
        
        fclose(fin_main);
        
    end % end rcut loop
