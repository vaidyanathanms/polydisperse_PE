%% To check equilibration
% Assumes chains have IDs from (nch_graft+1) to (nch_freearr+nch_graft)

clc
clear
close all
format long

%% Input data

%% Color codes
green = [0 0.5 0.0]; gold = [0.9 0.75 0]; orange = [0.91 0.41 0.17];brown = [0.2 0 0];
pclr = {'m',brown,green,'k','m', gold};
lsty = {'-','--',':'};
msty = {'d','s','o','x'};
cmap = [.7 .7 .7 %// light gray
    1  1  1]; %// white

%% Inputs
nch_graft = 32;
nch_freearr = [150];
nchain_mws = [52,8]; % Give the two molecular weights that correspond to the bidispersed case (larger MW first)
casearr  = [1];
pdi_freearr = [1.5];
arch_arr = {'bl_bl';'al_al'};
leg_arr  = {'Block-Block';'Alter-Alter'}; % ALWAYS CHECK for correspondence with arch_arr for legends
pdigraft = 1.0;
pdigraft_str = num2str(pdigraft,'%1.1f');
cutoff = '1.50';
ref_mw = 52;
dtval  = 12.5;

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
                
                tcnt = 1; tref = alldata.data(1,1); tarr_cnt = 1;
                for tcnt = 1:lendata
                    if alldata.data(tcnt,1) ~= tref
                        tarr_cnt = tarr_cnt + 1;
                        tref = alldata.data(tcnt,1);
                    end
                    if alldata.data(tcnt,5) == ref_mw
                        mol_id = alldata.data(tcnt,4) - nch_graft; %to make values
                        ads_matrix(mol_id,tarr_cnt) = 1;
                    end
                end
                
                yticks = 0:8:76;
                [r,c] = size(ads_matrix);
                
                h= figure;
                imagesc(dtval*((1:c)+0.5), (1:r)+0.5, ads_matrix) %https://stackoverflow.com/questions/3280705/how-can-i-display-a-2d-binary-matrix-as-a-black-white-plot
                colormap(cmap)
                ylim([1 75])
                xlim([0 (c+0.5)*dtval])
                axis ij
                set(gca,'YMinorGrid','On','YTick',yticks,'YTickLabel',yticks,'YGrid','On','FontSize',16)
                xlabel('Time ($\tau$)','FontSize',20,'Interpreter','Latex')
                ylabel('Chain ID','FontSize',20,'Interpreter','Latex')
                saveas(h,sprintf('./../../all_figures/bidispequildist_case_%d_arch_%s.png',casenum,dirstr));

            end % end case loop
            
        end % end arch loop
        
    end % end pdi loop
    
end % end nfree_arr loop
