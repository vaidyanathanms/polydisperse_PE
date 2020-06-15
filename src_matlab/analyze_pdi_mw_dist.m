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
nch_freearr = [16;32;64;96;128;150];
casearr  = [1;2;3;4];
pdi_freearr = [1.5];
arch_arr = {'bl_bl';'bl_al';'al_bl';'al_al'};
leg_arr  = {'Block-Block';'Block-Alter';'Alter-Block';'Alter-Alter'}; % ALWAYS CHECK for correspondence with arch_arr for legends
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
num_cases = length(casearr);

%% Main Analysis

if pdiflag % create consolidated list
    fout_cons = fopen(sprintf('./../../outfiles/overall/pdi_consolidate_rcut_%s.dat',...
        cutoff),'w');
    fprintf(fout_cons,'%s\t%s\t%s\t%s\t%s\t%s\n','N_f','Arch','Case_num',...
        'Mw_free','Mn_free','PDI_free');
end

for pdi_cntr = 1:length(pdi_freearr) % begin pdi free loop
    pdifree     = pdi_freearr(pdi_cntr);
    pdifree_str = num2str(pdifree,'%1.1f');
    
    %zero arrays for averages across cases
    casecntr_arr  = zeros(length(nch_freearr),length(arch_arr));
    npdi_all = zeros(length(nch_freearr),length(arch_arr));
    
    if pdiflag %create average across all cases
        favg_pdi = fopen(sprintf('./../../outfiles/overall/pdi_ave_allcases_rcut_%s.dat',...
            cutoff),'w');
        fprintf(favg_pdi,'%s\t%s\t%s\t%s\n','N_f','Arch','ncases','pdi_avg');
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
            
            
            if mwdflag %create average across all cases
                avgmolarr = zeros(num_cases*nval,1); % to compute average distribution
                favg_dist = fopen(sprintf('./../../outfiles/overall/avgdist_allcases_%s_rcut_%s.dat',...
                    dirstr,cutoff),'w');
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
                
                molarr = analyze_datafile(ads_fylename,nval,nch_graft); % extract molecular details
                
                if pdiflag % begin pdi calculation
                    
                    % analyze and write molecule details
                    moloutfyle = strcat(dirname,'/','init_mol_details.dat');
                    fmol = fopen(moloutfyle,'w');
                    fprintf(fmol,'%s\t%s\t%s\n','Actual_MolID','Remapped_MolID', 'Mol_Wt');
                    fprintf(fmol,'%d\t%d\t%d\n',[molarr(:,1) molarr(:,2) molarr(:,3)]');
                    fclose(fmol);
                    
                    % compute and write the initial pdi. Store avg arrays
                    [pdifree,mnfree,mwfree] = compute_pdi(molarr,nval);
                    fprintf(fout_cons,'%d\t%s\t%d\t%g\t%g\t%g\n',nval,dirstr,casenum,mnfree,mwfree,pdifree);
                    casecntr_arr(ncnt,arch_cnt)  = casecntr_arr(ncnt,arch_cnt) + 1;
                    npdi_all(ncnt,arch_cnt) = npdi_all(ncnt,arch_cnt) + pdifree;
                end % end pdi calculation
                
                if mwdflag % begin molecular weight distribution calculation
                    
                    %compute_mwdist();
                    outdist = compute_mwdist(molarr,3); % compute and write individual distribution
                    initdistoutfyle = strcat(dirname,'/','init_distout_details.dat');
                    fdist = fopen(initdistoutfyle,'w');
                    left_edges = outdist.BinEdges(1:outdist.NumBins);
                    right_edges = outdist.BinEdges(2:outdist.NumBins+1);
                    fprintf(fdist,'%s\t%d\n','NumBins',outdist.NumBins);
                    fprintf(fdist,'%s\t%s\t%s\n','leftEdge','RightEdge', 'Normalized Value');
                    fprintf(fdist,'%d\t%d\t%d\n',[left_edges' right_edges' outdist.Values']');
                    fclose(fdist);
                    
                    % Store into avg arrage
                    minmolval = 1+nval*(casecntr-1); maxmolval = casecntr*nval;
                    avgmolarr(minmolval:maxmolval,1) = molarr(:,3);
                    
                end % end mol. wt distribution calculation
                
                clear molarr
                
            end % end case loop
            
            if pdiflag
                % write avg pdi across cases
                avg_across_cases(ncnt,arch_cnt,pdi_cntr) = npdi_all(ncnt,arch_cnt)/casecntr_arr(ncnt,arch_cnt);
                fprintf(favg_pdi,'%d\t%s\t%d\t%g\n',nval,dirstr,casecntr_arr(ncnt,arch_cnt),avg_across_cases(ncnt,arch_cnt,pdi_cntr));
            end
            
            if mwdflag
                % find avg mw distribution across cases
                avgoutdist = compute_mwdist(avgmolarr,1);
                left_edges = avgoutdist.BinEdges(1:avgoutdist.NumBins);
                right_edges = avgoutdist.BinEdges(2:avgoutdist.NumBins+1);
                fprintf(favg_dist,'%s\t%d\n','NumBins',avgoutdist.NumBins);
                fprintf(favg_dist,'%s\t%s\t%s\n','leftEdge','RightEdge', 'Normalized Value');
                fprintf(favg_dist,'%d\t%d\t%d\n',[left_edges' right_edges' avgoutdist.Values']');
                fclose(favg_dist);
            end
            
        end % end arch loop
        
    end % end nfree loop
    
    fclose(favg_pdi);
    
end % end pdi free loop

fclose(fout_cons);