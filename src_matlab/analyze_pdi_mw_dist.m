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
nch_freearr = [16]%;32;64;96;128;150];
casearr  = [1]%,2,3,4];
pdi_freearr = [1.5];
arch_arr = {'bl_bl','bl_al','al_bl','al_al'};
leg_arr  = {'Block-Block','Block-Alter','Alter-Block','Alter-Alter'}; % ALWAYS CHECK for correspondence with arch_arr
pdigraft = 1.0;
cutoff = '1.50';
lz = 120; area=35^2;
set_tmax = 3e7; % maximum timestep for analysis;
nfreemons = 30; % this is average value or for PDI = 1.0 only.

%% Graft details

ncharge_mons_graft = 30; % per graft chain details
ntail_mons_graft = 5;    % per graft chain details
ntot_mons_graft  = ncharge_mons_graft + ntail_mons_graft; 
nch_graft = 32; 

%% Input flags
pdiflag  = 1;
mwdflag  = 1;

%% Zero arrays
avg_across_cases = zeros(length(nch_freearr),length(arch_arr),length(pdi_freearr));

%% Pre-calculations
rhofree = nch_freearr*nfreemons/(lz*area);
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
    casecntr_arr  = zeros(length(nch_freearr),length(arch_arr));
    nadschain_all = zeros(length(nch_freearr),length(arch_arr));
    totsamples    = zeros(length(nch_freearr),length(arch_arr));
    
    if pdiflag %create average across all cases
        fout_avg = fopen(sprintf('./../../outfiles/overall/pdi_ave_allcases_rcut_%s.dat',...
            cutoff),'w');
        fprintf(fout_avg,'%s\t%s\t%s\t%s\n','N_f','Arch','ncases','pdi_avg');
    end
    
    for ncnt = 1:length(nch_freearr) % begin nfree loop
        nval = nch_freearr(ncnt);
        
        for arch_cnt = 1:length(arch_arr)  % begin arch loop
            dirstr = arch_arr{arch_cnt};
            
            dirname = sprintf('./../../data_all_dir/n_%d/%s/pdifree%s_pdigraft_%s',...
                nval,dirstr,pdifree_str,pdigraft_str);
            
            if ~exist(dirname,'dir')
                fprintf('%s does not exist\n',dirname);
                continue
            end
            
            for casecntr = 1:length(casearr) % begin case loop
                casenum = casearr(casecntr);
                
                dirname = sprintf('./../../data_all_dir/n_%d/%s/pdifree%s_pdigraft_%s/Case_%d',...
                    nval,dirstr,pdifree_str,pdigraft_str,casenum);
                
                if ~exist(dirname,'dir')
                    fprintf('%s does not exist\n',dirname);
                    continue
                end
                
                fprintf('Analyzing pdi/nfree/arch/case: %g\t%d\t%s\t%d\n', pdifree,nval,dirstr,casenum);
                if pdiflag % begin adsorption calculation
                    
                    %check if file type exists
                    ads_prefix = sprintf('PEinitdata.txt');
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
                    
                    % begin running through the datafile
                    ads_fylename = strcat(dirname,'/',ads_fylelist(1).name);
                    if exist(ads_fylename,'file') ~= 2
                        fprintf('%s does not exist/empty file\n',ads_fylename);
                        continue;
                    elseif struct(dir(ads_fylename)).bytes == 0
                        fprintf('Empty file: %s \n',ads_fylename);
                        continue;
                    end
                    
                    % analyze and write molecule details
                    molarr_cnt = analyze_datafile(ads_fylename,nval,nch_graft);
                    moloutfyle = strcat(dirname,'/','init_mol_details.dat');
                    fmol = fopen(moloutfyle,'w');
                    fprintf(fmol,'%s\t%s\t%s\n','Actual_MolID','Remapped_MolID', 'Mol_Wt');
                    fprintf(fmol,'%d\t%d\t%d\n',[molarr_cnt(:,1) molarr_cnt(:,2) molarr_cnt(:,3)]');
                    fclose(fmol);
                    
                    
                    %compute_init_pdi();
                    %compute_mwdist();
                    
                    %save it to overall arrays
                    %casecntr_arr(ncnt,arch_cnt)  = casecntr_arr(ncnt,arch_cnt) + 1;
                    %nadschain_all(ncnt,arch_cnt) = nadschain_all(ncnt,arch_cnt) + avg_for_each_casenum;
                    %totsamples(ncnt,arch_cnt)    = totsamples(ncnt,arch_cnt) + tot_cntr_across_files;
                    
                end %end adsorption calculation
                
            end % end case loop
            
%             
%             avg_across_cases(ncnt,arch_cnt,pdi_cntr) = nadschain_all(ncnt,arch_cnt)/casecntr_arr(ncnt,arch_cnt);
%             fprintf(fout_avg,'%d\t%s\t%d\t%d\t%g\n',nval,dirstr,casecntr_arr(ncnt,arch_cnt),totsamples(ncnt,arch_cnt),avg_across_cases(ncnt,arch_cnt,pdi_cntr));
%             
        end % end arch loop
        
    end % end nfree loop
    
%     fclose(fout_avg);
    
end % end pdi free loop

fclose(fout_cons);