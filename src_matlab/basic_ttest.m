%% To compare between architectures and within the same architecture for correlation

% V2.0: Use ref1 and read like RDF calculations. For a fixed n_f and pdi_free and check for different
% correlations). Read from ttest_dir
% Added two sample unequal variance t-test (to be done: one sample)
% USE TAB AS DELIMITER WHEN USING TEXT TO COLUMNS IN EXCEL TO CHECK
% OUTPUTS. (If it is written as one column).

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
nfree_arr = [16;32;64;96;128;150];
max_numcases = 4; % how many MAXIMUM cases are available per n_pa
pdi_freearr = [1.5];
ref_arch_arr1 = {'bl_bl','bl_al','al_bl','al_al'};
pdigraft = 1.0;
nmonfree = 30; nmongraft = 30; ngraft = 32;
cutoff_arr = {'1.50'};
lz = 120; area=35^2;
set_tmax = 3e7; % maximum timestep for analysis;

%% Pre-calculations
rhofree = nfree_arr*30/(lz*area);
pdigraft_str = num2str(pdigraft,'%1.1f');
err_tol = 1e-10; % for finding elements in a real array

%% Main Analysis
for rcutcntr = 1:length(cutoff_arr) % begin rcut loop
    cutoff = cutoff_arr{rcutcntr};
    fprintf('Analyzing rcut = %s\n', cutoff);
    fylename = sprintf(sprintf('./../../ttest_dir/overall/alltruehypothesis_rcut_%s.dat',...
        cutoff));
    fout_true = fopen(fylename,'w');
    
    fprintf(fout_true,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n',...
        'nval','pdival','dirstr1','num_cases1','dirstr2','num_cases2','nullHyp','tvalue');
    
    for ncntr = 1:length(nfree_arr) % begin nfree loop
        nval = nfree_arr(ncntr);
        % write t-test file
        fylename = sprintf(sprintf('./../../ttest_dir/overall/ttestvals_n_%d_rcut_%s.dat',...
            nval,cutoff));
        fout_compare = fopen(fylename,'w');
        
        %out format:
        %pdival arch1 c1 .. c4 cavg tab tab arch2 tab tab arch2 c1 .. c4 cavg
        %tab tab tvalue
        
        fprintf(fout_compare,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t\t%s\t%s\t%s\t%s\t%s\t%s\t\t%s\t%s\t%s\t%s\n',...
            'nval','pdival','arch1','c1','c2','c3','c4','cavg',...
            'arch2','c1','c2','c3','c4','cavg','nullHyp','tvalue','ci_lo','ci_hi');
        
        for pdicntr = 1:length(pdi_freearr) % begin pdival loop
            pdival = pdi_freearr(pdicntr);
            
            
            for arch_cntr_1 = 1:length(ref_arch_arr1)-1 % begin ref_arch1 first loop
                
                % read file written by adsfrac for dirstr1 (arch_cntr_1)
                dirstr1 = ref_arch_arr1{arch_cntr_1};
                fylename = sprintf('./../../ttest_dir/n_%d/adsfrac_rcut_%s_pdifree_%g_arch_%s.dat',...
                    nval,cutoff,pdival,dirstr1);
                fin_main = fopen(fylename,'r'); % use same file ID so that no two files are opened at the same time to avoid confusion.
                
                if fin_main <= 0 % check for average list
                    fprintf('%s does not exist', fylename);
                    continue;
                end
                
                fgetl(fin_main); %skip first line
                tline = fgetl(fin_main);
                spl_tline = strsplit(strtrim(tline));
                len_tline1 = length(spl_tline); % use this as a separate number for writing into consolidated file
                avg_fvals_ref1 = zeros(len_tline1,1); % corresponding to num of cases in the file.
                
                for colcntr = 1:len_tline1
                    avg_fvals_ref1(colcntr)   = str2double(spl_tline{colcntr});
                end
                
                fclose(fin_main); % CLOSE fin_main
                clear tline spl_tline% clear this so that it is not overwritten
                %                 [hnull,pnull,cinull,statsnull] = ttest(avg_fvals_ref1); % ttest for equal variances and from same sample
                
                for arch_cntr_2 = 1+arch_cntr_1:length(ref_arch_arr1) % begin ref_arch1 second loop
            
                    % Do all the file out operations for the first
                    % architecture before moving to the second one.
                    fprintf(fout_compare,'%g\t',nval);
                    fprintf(fout_compare,'%g\t',pdival);
                    fprintf(fout_compare,'%s\t',dirstr1);
                    for colcntr = 1:len_tline1
                        fprintf(fout_compare,'%g\t',avg_fvals_ref1(colcntr)); % the values are already stored.
                    end
                    if len_tline1 ~= max_numcases % to fill the uncomputed cases with tabs so that it can be viewed easily in EXCEL sheets.
                        for fill_cntr = 1:max_numcases-len_tline1
                            fprintf(fout_compare,'\t');
                        end
                    end
                    fprintf(fout_compare,'%g\t \t', mean(avg_fvals_ref1)); % write mean value
                    
                    % Now read dirstr2 data (arch_cntr_2)
                    dirstr2 = ref_arch_arr1{arch_cntr_2};
                    fprintf('Analyzing n = %d, pdi = %g, arch1 = %s, arch2 = %s \n', ...
                        nval, pdival, dirstr1, dirstr2);
                    
                    % read file written by adsfrac for dirstr2
                    fylename = sprintf('./../../ttest_dir/n_%d/adsfrac_rcut_%s_pdifree_%g_arch_%s.dat',...
                        nval,cutoff,pdival,dirstr2);
                    fin_main = fopen(fylename,'r');
                    
                    if fin_main <= 0 % check for average list
                        fprintf('%s does not exist', fylename);
                        continue;
                    end
                    
                    fgetl(fin_main); %skip first line
                    tline = fgetl(fin_main);
                    spl_tline = strsplit(strtrim(tline));
                    len_tline2 = length(spl_tline); % should be different so that len_tline1 is not overwritten for outfile requirements
                    avg_fvals_ref2 = zeros(len_tline2,1); % corresponding to num of cases in the file.
                    
                    fprintf(fout_compare,'%s\t',dirstr2);
                    for colcntr = 1:len_tline2
                        avg_fvals_ref2(colcntr) = str2double(spl_tline{colcntr});
                        fprintf(fout_compare,'%g\t',avg_fvals_ref2(colcntr));
                    end
                    
                    if len_tline2 ~= max_numcases % to fill the uncomputed cases with tabs so that it can be viewed easily in EXCEL sheets.
                        for fill_cntr = 1:max_numcases-len_tline2
                            fprintf(fout_compare,'\t');
                        end
                    end
                    fprintf(fout_compare,'%g\t \t', mean(avg_fvals_ref2)); % write mean value
                    fclose(fin_main); % CLOSE fin_main
                    clear tline spl_tline len_tline2% clear this so that it is not overwritten
                    
                    [hnull,pnull,cinull,statsnull] = ttest2(avg_fvals_ref1,avg_fvals_ref2,'Vartype','unequal'); % ttest for unequal variances and different samples
                    
                    fprintf(fout_compare,'%d\t%g\t%g\t%g\n',hnull,pnull,cinull(1),cinull(2)); % write the ttest results
                    
                    if hnull
                        fprintf(fout_true,'%d\t%g\t%s\t%d\t%g\t%s\t%d\t%g\t%d\t%g\n',...
                            nval,pdival,dirstr1,len_tline1,dirstr2,len_tline2,hnull,pnull);
                    end
                    
                end % end arch_cntr_2 (dirstr2)
                
                clear len_tline1
                
            end % end arch_cntr_1 (dirstr1)
   
        end % end pdi_cntr
            
        fclose(fout_compare); % close outfile

    end % end nval loop

    fclose(fout_true);
    
end % end rcut loop
