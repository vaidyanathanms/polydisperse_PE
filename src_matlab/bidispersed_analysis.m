%% To analyze the bidispersed cases

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

%% Flags
plt_data = 1; % Plot ratio of histogram

%% Zero arrays

mean_data_arch_rat   = zeros(length(arch_arr));
mean_data_arch_big   = zeros(length(arch_arr));
mean_data_arch_sml   = zeros(length(arch_arr));

stderr_data_arch_rat = zeros(length(arch_arr));
stderr_data_arch_big = zeros(length(arch_arr));
stderr_data_arch_sml = zeros(length(arch_arr));

%% Main Analysis

for ncnt = 1:length(nch_freearr) % begin nfree loop
    nval = nch_freearr(ncnt);
    
    % create consolidated list
    fout_cons = fopen(sprintf('./../../bidispersed_cases/fracads_n_%d_rcut_%s.dat',...
        nval,cutoff),'w');
    fprintf(fout_cons,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n','N_f','Arch','Ref_PDI','Case_num',...
        'nsteps','Mw1_free','Mw2_free','fads_cnt_mw1','fads_cnt_mw2','rat_cnt_fads');
    
    %create average across all cases
    favg_pdi = fopen(sprintf('./../../bidispersed_cases/fracads_ave_allcases_rcut_%s.dat',...
        cutoff),'w');
    fprintf(favg_pdi,'%s\t%s\t%s\t%s\t%s\t%s\n','N_f','Arch','PDI_free','ncases','mean_rat_fads','SEM');
    
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
                
                
                mwbigarr(casecntr,1) = 2*mwcnt1/(tsteps*nval); % ratio wrt initial total steps * number of chains of each type (half of each type and hence 2)
                mwsmlarr(casecntr,1) = 2*mwcnt2/(tsteps*nval); % ratio wrt initial total steps * number of chains of each type (half of each type and hence 2)
                
                fprintf(fout_cons,'%d\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%g\n', nval,dirstr,pdifree_str,casenum,...
                    tsteps,nchain_mws(1),nchain_mws(2),mwcnt1,mwcnt2,mwrat);
                
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
            
            
            mean_data_arch_rat(arch_cnt,1)   = meanratval;
            mean_data_arch_big(arch_cnt,1)   = meanbigval;
            mean_data_arch_sml(arch_cnt,1)   = meansmlval;
            
            stderr_data_arch_rat(arch_cnt,1) = stdraterr;
            stderr_data_arch_big(arch_cnt,1) = stdbigerr;
            stderr_data_arch_sml(arch_cnt,1) = stdsmlerr;
            
        end % end arch loop
        
    end % end pdi loop
    
    fclose(fout_cons);
    fclose(favg_pdi);
    
end % end nfree_arr loop


%% Plot histogram


if plt_data
    
    % Plot only ratio
    h1 = figure;
    hold on
    box on
    set(gca,'FontSize',16)
    xlabel('Architecture','FontSize',20,'Interpreter','Latex')
    ylabel('$\phi_{\rm{ads}}$','FontSize',20,'Interpreter','Latex')
    X = categorical({leg_arr{1},leg_arr{2}});
    bar(X,mean_data_arch_rat(:,1));
    er = errorbar(X,mean_data_arch_rat(:,1),stderr_data_arch_rat(:,1),stderr_data_arch_rat(:,1));
    er.Color = [0 0 0];
    er.LineStyle = 'none';
    er.LineWidth = 2;
    saveas(h1,'./../../Figs_paper/bidispersed.png');
    
    % Plot all mws and ratio of mws
    h2 = figure;
    hold on
    box on
    set(gca,'FontSize',16)
    xlabel('Architecture','FontSize',20,'Interpreter','Latex')
    ylabel('$f_{\rm{ads}}$','FontSize',20,'Interpreter','Latex')
    X = categorical({leg_arr{1};leg_arr{2}});
    Yall_data = zeros(length(arch_arr),3); Erall_data = zeros(length(arch_arr),3);
    for i = 1:length(arch_arr)
        Yall_data(i,1) = mean_data_arch_sml(i,1);
        Yall_data(i,2) = mean_data_arch_big(i,1);
        Yall_data(i,3) = mean_data_arch_rat(i,1);
        
        Erall_data(i,1) = stderr_data_arch_sml(i,1);
        Erall_data(i,2) = stderr_data_arch_big(i,1);
        Erall_data(i,3) = stderr_data_arch_rat(i,1);
    end

    
    hBar = bar(X,Yall_data);
    ylim([0.006 15])
    pleg{1} = '$n$ = 8';
    pleg{2} = '$n$ = 52';
    pleg{3} = '$\phi_{\rm{ads}}$';
    legend(hBar,pleg,'Interpreter','Latex','Location','NorthWest')
    %xticks([1 2]); xticklabels = ({leg_arr{1} leg_arr{2}});
    
    % Only to plot with error bars
    % Ref: https://www.mathworks.com/matlabcentral/answers/102220-how-do-i-place-errorbars-on-my-grouped-bar-graph-using-function-errorbar-in-matlab
%     % Find the number of groups and the number of bars in each group
%     ngroups = size(Yall_data, 1);
%     nbars = size(Yall_data, 2);
%     % Calculate the width for each bar group
%     groupwidth = min(0.8, nbars/(nbars + 1.5));
%     % Set the position of each error bar in the centre of the main bar
%     % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
%     for i = 1:nbars
%         % Calculate center of each bar
%         x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
%         errorbar(x, Yall_data(:,i), Erall_data(:,i), 'k', 'linestyle', 'none');
%     end
    set(gca,'yscale','log'); % edit legends and xlabels AFTER PLOTTING. Could not insert automatically.
    hold off
    
   
end


