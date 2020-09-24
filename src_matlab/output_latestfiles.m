%% To output the latest files along with whether the analysis has been finished
% Dependencies: find_latest_trfyle.m; find_latest_anafyle.m

clc;
clear;
close all;
format long

%% Color codes
green = [0 0.5 0.0]; gold = [0.9 0.75 0]; orange = [0.91 0.41 0.17];brown = [0.2 0 0];
pclr = {'m',brown,green,'k','m', gold};
lsty = {'-','--',':'};
msty = {'d','s','o','x'};

%% Inputs
nch_freearr = [16]%;32;64;128];
casearr  = [1]%;2;3;4];
pdi_freearr = [1.5];
arch_arr = {'bl_bl';'al_al'};
leg_arr  = {'Block-Block';'Alter-Alter'}; % ALWAYS CHECK for correspondence with arch_arr for legends
pdigraft = 1.0;
pdigraft_str = num2str(pdigraft,'%1.1f');
cutoff = '1.50';

%% Begin analysis

for pdi_cntr = 1:length(pdi_freearr) % begin pdi free loop
    
    ref_pdifree     = pdi_freearr(pdi_cntr);
    pdifree_str     = num2str(ref_pdifree,'%1.1f');

    dirname = sprintf('./../../all_trajfiles');
    if ~exist(dirname,'dir')
        fprintf('%s does not exist\n',dirname);
        continue
    end

    fout_cons = fopen(sprintf('./../../all_trajfiles/latest_pdi_%s.dat',pdifree_str),'w');    
    fprintf(fout_cons,'%s\t%s\t%s\t%s\t%s\t%s\t%s\n','N_pa','Arch','Casenum','latest_file',...
        'denscalc','adscalc','mwdistcalc');
    
    % Read trajlist
    trlist_fname = strcat(dirname,'/',sprintf('all_trajfiles_pdi_%s.txt',pdifree_str));
    if exist(trlist_fname,'file') ~= 2
        fprintf('%s does not exist/empty file\n',trlist_fname);
        continue;
    else
        flist = fopen(trlist_fname);
        tline = fgetl(flist); %skip first line
        lcnt = 0;
        alldata = cell(1,5);
        while ~feof(flist)
            tline = fgetl(flist);
            if ~ischar(tline) || isempty(tline)
                continue;
            end
            spl_tline = strsplit(strtrim(tline));
            if length(spl_tline) ~= 5
                fprintf('ERR: Not all data present in the line of %s: %s',trlist_fname,spl_tline);
                continue;
            end
            lcnt = lcnt + 1;
            for wcnt = 1:5
                alldata{lcnt,wcnt} = spl_tline{wcnt};
            end
        end
    end
    
    for ncnt = 1:length(nch_freearr) % begin nfree loop
        nval = nch_freearr(ncnt);
    
        for arch_cnt = 1:length(arch_arr)  % begin arch loop
            dirstr = arch_arr{arch_cnt};
            
            for casecntr = 1:length(casearr) % begin case loop
                casenum = casearr(casecntr);
                traj_fylelist = cell(1,1);
                fylcnt = 0;
                
                fprintf('Analyzing pdi/Npa/arch/casenum %s\t%d\t%s\t%d\n',pdifree_str,nval,dirstr,casenum);
                
                for fcnt = 1:length(alldata(:,1))

                    % Check for trajfile
                    if nval == str2double(alldata(fcnt,1)) && strcmp(dirstr,alldata(fcnt,2)) ...
                            && casenum == str2double(alldata(fcnt,4))
                        fylcnt = fylcnt + 1;
                        traj_fylelist{fylcnt,1} = alldata(fcnt,5);
                    end
                    
                end % end iterating through all files for the given condition
                
                if min(size(traj_fylelist)) == 0 || numel(traj_fylelist) == 0 || isempty(traj_fylelist{1,1})
                    fprintf('No trajfiles are found for %d\t%s\t%d\n',nval,dirstr,casenum);
                    fprintf(fout_cons,'%d\t%s\t%d\t%s\n',nval,dirstr,casenum,'UNK');
                else
                    latest_index = find_latest_trfyle(traj_fylelist,length(traj_fylelist),2);
                    if latest_index == -1
                        fprintf(fout_cons,'%d\t%s\t%d\t%s\t',nval,dirstr,casenum,'UNK');
                    else
                        fprintf(fout_cons,'%d\t%s\t%d\t%s\t',nval,dirstr,casenum,char(traj_fylelist{latest_index}));
                    end
                end
                
                % Check for grpdens
                ana_dirname = sprintf('./../../sim_results/outresults_dir_n_%d/%s/pdifree_%s_pdigraft_%s/Case_%d',...
                    nval,dirstr,pdifree_str,pdigraft_str,casenum);
                if ~exist(ana_dirname,'dir')
                    fprintf('%s does not exist\n',ana_dirname);
                    continue
                end
                dens_prefix = 'dens_config_*.lammpstrj';
                dens_fylelist = dir(strcat(ana_dirname,'/',dens_prefix));
                if min(size(dens_fylelist)) == 0 || numel(dens_fylelist) == 0
                    fprintf('No files/Empty files are found for %s\n',dens_prefix);
                    fprintf(fout_cons,'%s\t','N/A');
                else
                    [latest_index] = find_latest_anafyle(dens_fylelist,numel(dens_fylelist),3);
                    if latest_index == -1
                        fprintf(fout_cons,'%s\t','UNK');
                    else
                        fprintf(fout_cons,'%s\t',dens_fylelist(latest_index).name);
                    end
                end
                
                % Check for adsfrac
                ads_prefix = sprintf('adsfracchain_rcut_%s_config_*.lammpstrj',cutoff);
                ads_fylelist = dir(strcat(ana_dirname,'/',ads_prefix));
                if min(size(ads_fylelist)) == 0 || numel(ads_fylelist) == 0
                    fprintf('No files/Empty files are found for %s\n',ads_prefix);
                    fprintf(fout_cons,'%s\t','N/A');
                else
                    [latest_index] = find_latest_anafyle(ads_fylelist,numel(ads_fylelist),6);
                    if latest_index == -1
                        fprintf(fout_cons,'%s\t','UNK');
                    else
                        fprintf(fout_cons,'%s\t',ads_fylelist(latest_index).name);
                    end
                end
                
                %check for mwdist
                mw_prefix = sprintf('chainsadsval_rcut_%s_config_*.lammpstrj',cutoff);
                mw_fylelist = dir(strcat(ana_dirname,'/',mw_prefix));
                if min(size(mw_fylelist)) == 0 || numel(mw_fylelist) == 0
                    fprintf('No files/Empty files are found for %s\n',mw_prefix);
                    fprintf(fout_cons,'%s\t','N/A');
                else
                    [latest_index] = find_latest_anafyle(mw_fylelist,numel(mw_fylelist),6);
                    if latest_index == -1
                        fprintf(fout_cons,'%s\t','UNK');
                    else
                        fprintf(fout_cons,'%s\t',mw_fylelist(latest_index).name);
                    end
                end
                
            end % end casenum
            
        end % end arch_cnt
        
    end % end ncnt
    
end % end pdifree

                
            
    