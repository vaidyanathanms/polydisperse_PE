%% To convert dens*/grpdens* to polydens*
%% This is a workaround.
% This code assumes that the density profiles of the free chains are
% identical to that in grpdens*. For the graftchains, it can be found using
% \rho_g(z) = int(\rho_A(z)*NA+\rho_B(z)*NB)/(NA+NB) where NA and NB are
% neutral and charged monomers of the graft which is same for all the cases
% in this work.

clc;
clear;
close all;
format long;

%% Color codes
green = [0 0.5 0.0]; gold = [0.9 0.75 0]; orange = [0.91 0.41 0.17];brown = [0.2 0 0];
pclr = {'r','b',green,brown,'k','m', gold};
lsty = {'--','-',':'};
msty = {'d','o','s','x'};

%% Input data
nplot = 32;
casearr  = [1]%,3,4];
pdi_freearr = [1.5];
arch_arr = {'bl_bl','al_al'};
leg_arr  = {'Block-Block','Alter-Alter'}; % ALWAYS CHECK for correspondence with arch_arr
overwrite = 1; %BE VERY CAREFUL. THIS WILL OVERWRITE. USE WITH CAUTION
pdigraft = 1.0;
nmonfree = 30; nmongraft = 30; ngraft = 32;
cutoff = '1.50';
nch_graft = 32;
lz = 120; area=35^2;
set_tmax = 3e7; % maximum timestep for analysis;
tot_graftmon = nmongraft*ngraft;


%% Pre-calculations
pdigraft_str = num2str(pdigraft,'%1.1f');
fprintf('%s\n','Preparing density plots');

for arch_cnt = 1:length(arch_arr)  % begin arch loop
    dirstr = arch_arr{arch_cnt};lcnt = 1;
    
    % Plot density profiles
    h1 = figure;
    hold on
    box on
    set(gca,'FontSize',16)
    xlabel('$z/L_{z}$','FontSize',20,'Interpreter','Latex')
    ylabel('$\rho$','FontSize',20,'Interpreter','Latex')
    
    for pdi_cntr = 1:length(pdi_freearr) % begin pdi free loop
        
        pdifree     = pdi_freearr(pdi_cntr);
        pdifree_str = num2str(pdifree,'%1.1f');
        totcases = 0;avg_rho_graft = 0; avg_rho_free = 0;
        
        for casecntr = 1:length(casearr) % begin case loop
            casenum = casearr(casecntr);
            
            % Check file existence
            dirname = sprintf('./../../sim_results/outresults_dir_n_%d/%s/pdifree_%s_pdigraft_%s/Case_%d',...
                nplot,dirstr,pdifree_str,pdigraft_str,casenum);
            if ~exist(dirname,'dir')
                fprintf('%s does not exist\n',dirname);
                continue
            end
            
            % If file exists, DO NOT do anything
            rho_prefix = 'polydens_config_*.lammpstrj';
            rho_fylelist = dir(strcat(dirname,'/',rho_prefix));            
            nfyles = numel(rho_fylelist); %number of files of the type
            if nfyles ~= 0
                fprintf('Found files of the type polydens_config* in %d\t%s\n', nplot,dirstr);
                fprintf('WARNING: Better to use polydens_config* with plot_paper.m for plotting densities..\n')
            
                m = input('Do you want to continue, Y/N:\n','s');
                
                if m == 'N' || m == 'n'
                    continue;
                elseif m == 'Y' || m == 'y'
                    fprintf('CONTINUING. WILL OVERWRITE FILES\n')
                else
                    fprintf('ENTER EITHER Y or N\n')
                    error('CHECK INPUT\n')
                end
            end
            
            % Analyze grpdens_config_* to compute the free chain profiles    
            rho_prefix = 'grpdens_config_*.lammpstrj';
            rho_fylelist = dir(strcat(dirname,'/',rho_prefix));
            if min(size(rho_fylelist)) == 0
                fprintf('No files/Empty files are found for %s\n',rho_prefix);
                continue;
            end
                
            [latest_fyleindex] = find_latest_fyle(rho_fylelist,nfyles);
            fylename   = strsplit(rho_fylelist(latest_fyleindex).name,{'_','.'});
            tstamp     = str2double(fylename{3});
            freeindex  = tstamp;
             
            rho_fylename = strcat(dirname,'/',rho_fylelist(latest_fyleindex).name);
            if exist(rho_fylename,'file') ~= 2
                fprintf('%s does not exist/empty file\n',rho_fylename);
                continue;
            elseif struct(dir(rho_fylename)).bytes == 0
                fprintf('Empty file: %s \n',rho_fylename);
                continue;
            end
            all_data = importdata(rho_fylename,' ',1);
            fld = all_data.data;
            rdata = fld(:,1); 
            free_data = fld(:,3);
            nbin_free = length(rdata(:,1));
            
            % Analyze dens_* and covert to polydens using the formula given
            % at the top/header of this file.
            rho_prefix = 'dens_config_*.lammpstrj';
            rho_fylelist = dir(strcat(dirname,'/',rho_prefix));
            if min(size(rho_fylelist)) == 0
                fprintf('No files/Empty files are found for %s\n',rho_prefix);
                continue;
            end
                
            [latest_fyleindex] = find_latest_fyle(rho_fylelist,nfyles);
            fylename   = strsplit(rho_fylelist(latest_fyleindex).name,{'_','.'});
            tstamp     = str2double(fylename{3});
            graftindex = tstamp;
            
             
            rho_fylename = strcat(dirname,'/',rho_fylelist(latest_fyleindex).name);
            if exist(rho_fylename,'file') ~= 2
                fprintf('%s does not exist/empty file\n',rho_fylename);
                continue;
            elseif struct(dir(rho_fylename)).bytes == 0
                fprintf('Empty file: %s \n',rho_fylename);
                continue;
            end
            all_data = importdata(rho_fylename,' ',1);
            fld = all_data.data;
            rdata = fld(:,1);
            gr_neutral = fld(:,3); gr_charged = fld(:,4);
            nbin_graft = length(rdata(:,1));
            gr_all = 0.5*(gr_neutral+gr_charged); % same as 0.5*(gr_neutral*nmongraft+gr_charged*nmongraft)/nmongraft;
            
            % sanity check for length of data
            if nbin_free ~= nbin_graft
                fprintf('WARNING: The number of bins do not match between graft and free\n')
                error('CHECK nbins')
            elseif graftindex ~= freeindex
                fprintf('WARNING: Analyzing different time frames between graft and free\n')
                fprintf('WARNING: Using the graft time as the reference time for output file name\n')
            end
            
            
            out_prefix = strcat('polydens_config_',num2str(graftindex),'.lammpstrj');
            fout_file  = strcat(dirname,'/',out_prefix);
            fout_poly  = fopen(fout_file,'w');
            fprintf(fout_file,'w');
            for i = 1:length(rdata)
                fprintf(fout_poly,'%g\t%g\t%g\n',rdata(i,1),gr_all(i,1),free_data(i,1));
            end
            
        end
                
    end
    
end

