%% To check equilibration
% Assumes chains have IDs from 1-nch_freearr

clear
clc
close all
format long

%% Input data

%% Color codes
green = [0 0.5 0.0]; gold = [0.9 0.75 0]; orange = [0.91 0.41 0.17];brown = [0.2 0 0];
pclr = {'m',brown,green,'k','m', gold};
lsty = {'-','--',':'};
msty = {'d','s','o','x'};

%% Inputs
nch_freearr = [150];
nchain_mws = [52,8]; % Give the two molecular weights that correspond to the bidispersed case (larger MW first)
casearr  = [1;2;3;4];
pdi_freearr = [1.5];
arch_arr = {'bl_bl';'al_al'};
leg_arr  = {'Block-Block';'Alter-Alter'}; % ALWAYS CHECK for correspondence with arch_arr for legends
pdigraft = 1.0;
pdigraft_str = num2str(pdigraft,'%1.1f');
cutoff = '1.50';
ref_mw = 52;


%% Flags
plt_data = 1; % Plot ratio of histogram

%% Main Analysis

for ncnt = 1:length(nch_freearr) % begin nfree loop
    nval = nch_freearr(ncnt);
    
    for pdi_cntr = 1:length(pdi_freearr) % begin pdi free loop
        ref_pdifree     = pdi_freearr(pdi_cntr);
        pdifree_str     = num2str(ref_pdifree,'%1.1f');
        
        %zero arrays for averages across cases
        casecntr_arr  = zeros(length(nch_freearr),length(arch_arr));
        npdi_all = zeros(length(nch_freearr),length(arch_arr));
        
        for arch_cnt = 1:length(arch_arr)  % begin arch loop
            dirstr = arch_arr{arch_cnt};
            
            dirname = sprintf('./../../bidispersed_cases/outresults_dir_n_%d/%s/pdifree_%s_pdigraft_%s',...
                nval,dirstr,pdifree_str,pdigraft_str);
            
            if ~exist(dirname,'dir')
                fprintf('%s does not exist\n',dirname);
                continue
            end
            
            mwratarr = zeros(length(casearr),1);
            mwbigarr = zeros(length(casearr),1);
            mwsmlarr = zeros(length(casearr),1);
            
            for casecntr = 1:length(casearr) % begin case loop
                casenum = casearr(casecntr);
                
                dirname = sprintf('./../../bidispersed_cases/outresults_dir_n_%d/%s/pdifree_%s_pdigraft_%s/Case_%d',...
                    nval,dirstr,pdifree_str,pdigraft_str,casenum);
                
                if ~exist(dirname,'dir')
                    fprintf('%s does not exist\n',dirname);
                    continue
                end
                
                fprintf('Analyzing pdi/nfree/arch/case: %g\t%d\t%s\t%d\n', ref_pdifree,nval,dirstr,casenum);
                %check if file type exists
                ads_prefix = sprintf('chainadsval_rcut_%s_config_*.lammpstrj',cutoff);
                ads_fylelist = dir(strcat(dirname,'/',ads_prefix));
                if min(size(ads_fylelist)) == 0
                    fprintf('No files/Empty files are found for %s\n',ads_prefix);
                    continue;
                end
                
                len_adsfyles = numel(ads_fylelist); %number of files of the type
                ads_fylename = strcat(dirname,'/',ads_fylelist(len_adsfyles).name); % select the latest file
                
                fprintf('File under processing: %s\n', ads_fylename);
                
                alldata = importdata(ads_fylename);
                adsdata = alldata.data(:,5);
                lendata = length(adsdata(:,1));
                
                mwcnt1 = sum(adsdata(:) == nchain_mws(1));
                mwcnt2 = sum(adsdata(:) == nchain_mws(2));
                if (mwcnt1+mwcnt2) ~= lendata
                    fprintf('Error in cross checking total adsorbed molecules. Sum does not match \n');
                    fprintf('%d\t%d\t%d\n',mwcnt1,mwcnt2,lendata)
                    continue;
                end
                mwrat  = mwcnt1/mwcnt2;
                
                mwratarr(casecntr,1) = mwrat;
                
                % Count total time steps
                tsteps = 1; tref = alldata.data(1,1);
                for tcnt = 2:lendata
                    if tref ~= alldata.data(tcnt,1)
                        tsteps = tsteps + 1;
                        tref = alldata.data(tcnt,1);
                    end
                end
                
                % see whether chains are adsorbed/desorbed - Heart of the
                % code
                
                ads_matrix = zeros(nval, tsteps);
                
                mwbigarr(casecntr,1) = 2*mwcnt1/(tsteps*nval); % ratio wrt initial total steps * number of chains of each type (half of each type and hence 2)
                mwsmlarr(casecntr,1) = 2*mwcnt2/(tsteps*nval); % ratio wrt initial total steps * number of chains of each type (half of each type and hence 2)
                
                for molid_cnt = 1:nval
                    if molid_cnt
                    for tcnt = 1:tsteps
                        
                
                
            end
            
            % weed out zero values
            allnmwrat = mwratarr(mwratarr~=0);
            allnbigmw = mwbigarr(mwbigarr~=0);
            allnsmlmw = mwsmlarr(mwsmlarr~=0);
            
            lenarr = length(allnmwrat);
            stdratval = std(allnmwrat);
            meanratval = mean(allnmwrat);
            stdraterr = stdratval/sqrt(lenarr);
            
            fprintf(favg_pdi,'%d\t%s\t%s\t%d\t%g\t%g\n',nval,dirstr,pdifree_str,lenarr,meanratval,stdraterr);
            
            stdbigval = std(allnbigmw);
            meanbigval = mean(allnbigmw);
            stdbigerr = stdbigval/sqrt(lenarr);
            
            stdsmlval = std(allnsmlmw);
            meansmlval = mean(allnsmlmw);
            stdsmlerr = stdsmlval/sqrt(lenarr);
            

            
        end % end arch loop
        
    end % end pdi loop
    
    fclose(fout_cons);
    fclose(favg_pdi);
    
end % end nfree_arr loop
