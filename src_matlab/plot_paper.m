%% Plot all data for the paper

% Figsz   : SZ-distribution
% Figfads : f_ads vs Npa/Npc for different PDI and architectures
% Figqnet : Qnet vs Npa/Npc for different PDI and architectures
% Figdens : Density plots

clear
clc
close all
format long

%% Color codes
green = [0 0.5 0.0]; gold = [0.9 0.75 0]; orange = [0.91 0.41 0.17];brown = [0.2 0 0];
pclr = {'r','b',green,brown,'k','m', gold};
lsty = {'--','-',':'};
msty = {'d','o','s','x'};

%% Input data
nfreearr = [16;32;64;128;150];
casearr  = [1,2,3,4];
pdi_freearr = [1,1.5];
arch_arr = {'bl_bl','al_al'};
leg_arr  = {'Block-Block','Alter-Alter'}; % ALWAYS CHECK for correspondence with arch_arr
pdigraft = 1.0;
nmonfree = 30; nmongraft = 30; ngraft = 32;
cutoff = '1.50';
nch_graft = 32;
lz = 120; area=35^2;
set_tmax = 3e7; % maximum timestep for analysis;

%% Input flags
% see definitions above
figsz   = 0;
figfads = 1;
figqnet = 1;
figdens  = 0; nplot = 150; %for density profiles

%% Pre-calculations
rhofree = nfreearr*30/(lz*area);
pdigraft_str = num2str(pdigraft,'%1.1f');

%% SZ
if figsz
    
    nval_pl = [64,150]; casenum = 1; arch = 'bl_bl'; mnideal = 30; pdifree = 1.50;
    h1 = figure;
    hold on
    box on
    set(gca,'FontSize',16)
    xlabel('$M_{i}$','FontSize',20,'Interpreter','Latex')
    ylabel('$P(y)$','FontSize',20,'Interpreter','Latex')
    for plcnt = 1:length(nval_pl)
        dirname = sprintf('./../../data_all_dir/n_%d/%s/pdifree1.5_pdigraft_1.0/Case_%d',...
            nval_pl(plcnt),arch,casenum);
        alldata = importdata(strcat(dirname,'/init_mol_details.dat'));
        histogram(alldata.data(:,3),'BinWidth',8,'Normalization','pdf','BinLimits',[1,max(alldata.data(:,3))+1]);
        legendinfo{plcnt} = ['$n_{pa} =$ ' num2str(nval_pl(plcnt))];
    end
    
    mwvals = 1:10*mnideal;
    mwrat  = mwvals/mnideal;
    k = 1/(pdifree - 1);
    
    term1 = k^k;
    term2 = mwrat.^(k-1);
    term3 = exp(-k*mwrat);
    term4 = gamma(mwrat);
    
    psztheory = (term1.*term2.*term3)./term4;
    normval = trapz(mwvals,psztheory);
    
    
    plot(mwvals,psztheory/normval,'-k','LineWidth',2);
    legendinfo{plcnt+1} = 'Theory';
    xlim([0 max(alldata.data(:,3))+10])
    
    legend(legendinfo,'FontSize',16,'Location','Best','Interpreter','Latex')
    legend boxoff
    saveas(h1,'./../../Figs_paper/SZdistribution.png');
    clear legendinfo
    
end

%% f_ads
if figfads
    % Create arrays for according to the input array sizes
    frac_ads = zeros(length(nfreearr),length(arch_arr)+1,length(pdi_freearr)); %+1 for y-dimension to incorporate the N_f values
    err_ads  = zeros(length(nfreearr),length(arch_arr)+1,length(pdi_freearr)); %+1 for y-dimension to incorporate the N_f values
    pdi_plot = zeros(length(pdi_freearr));
    
    for pdi_cntr = 1:length(pdi_freearr) % begin pdi free loop
        
        pdifree     = pdi_freearr(pdi_cntr);
        pdifree_str = num2str(pdifree,'%1.1f');
        
        frac_ads(:,1,pdi_cntr) = nfreearr;
        err_ads(:,1,pdi_cntr)  = nfreearr;
        pdi_plot(pdi_cntr,1)   = pdifree;
        
        % Check file existence
        fname = sprintf('./../../outfiles/overall/adsorbed_chain_ave_allcases_rcut_%s_pdifree_%g.dat',...
            cutoff,pdifree);
        fads_id = fopen(fname);
        
        if fads_id <= 0
            fprintf('%s does not exist', fname);
            continue;
        end
        
        fprintf('File under process %s\n', fname);
        
        % Read and parse first line
        tline = fgetl(fads_id); % get header
        if ~ischar(tline) || isempty(tline)
            fprintf('ERROR: Unable to read file %s\n', tline);
            continue;
        end
        spl_tline = strtrim(strsplit(strtrim(tline)));
        
        nf_col = -1; arch_col = -1; favg_col = -1; err_col = -1;
        for wcnt = 1:length(spl_tline)
            if strcmp(spl_tline{wcnt},'N_f')
                nf_col = wcnt;
            elseif strcmp(spl_tline{wcnt},'Arch')
                arch_col = wcnt;
            elseif strcmp(spl_tline{wcnt},'avg_fraction')
                favg_col = wcnt;
            elseif strcmp(spl_tline{wcnt},'StdErrMean')
                err_col = wcnt;
            end
        end
        
        if nf_col == -1 || arch_col == -1 || favg_col == -1 || err_col == -1
            fprintf('Could not find all headers in %s \n' , tline);
            continue
        end
        
        
        while ~feof(fads_id)
            
            % Store data into arrays by comparing with the input arrays
            tline = fgetl(fads_id); % get line
            if ~ischar(tline) && isempty(tline)
                fprintf('ERROR: Unable to read line %s\n', tline);
                continue;
            end
            spl_tline = strtrim(strsplit(strtrim(tline)));
            
            % nf value
            findnf = -1;
            for itercnt = 1:length(nfreearr) % comparing wrt the INPUT array
                if str2double(spl_tline{nf_col}) == nfreearr(itercnt)
                    findnf = 1;
                    rownum = itercnt;
                    break;
                end
            end
            if findnf == -1
                fprintf('WARNING: Did not find N_pa value in the input array: %d\n', str2double(spl_tline{nf_col}));
                continue;
            end
            
            
            %arch value
            findarch = -1;
            for itercnt = 1:length(arch_arr) %comparing wrt the INPUT Array
                if strcmp(arch_arr{itercnt},spl_tline{arch_col})
                    findarch = 1;
                    colnum = itercnt;
                    break;
                end
            end
            
            if findarch == -1
                fprintf('WARNING: Did not find arch value in the input array: %s\n', spl_tline{arch_col});
                continue;
            end
            
            frac_ads(rownum,colnum+1,pdi_cntr) = str2double(spl_tline{favg_col});
            err_ads(rownum,colnum+1,pdi_cntr)  = str2double(spl_tline{err_col});
            
        end
        
        fprintf('Finished analyzing %s\n', fname);
        
    end
    
    % Plot frac_ads, err_ads in one go
    
    h1 = figure;
    hold on
    box on
    set(gca,'FontSize',16)
    xlabel('$n_{pa}/n_{pc}$','FontSize',20,'Interpreter','Latex')
    ylabel('$f_{\rm{ads}}$','FontSize',20,'Interpreter','Latex')
    
    lcnt = 1;
    for plcnt = 1:length(pdi_freearr)
        
        for arch_cnt = 1:length(arch_arr)
            
            if pdi_freearr(plcnt) == 1
                errorbar(frac_ads(:,1,plcnt)/nch_graft,frac_ads(:,1+arch_cnt,plcnt),err_ads(:,1+arch_cnt,plcnt),...
                    'Color', pclr{arch_cnt},'Marker',msty{arch_cnt},'MarkerFaceColor',...
                    'None','LineStyle',lsty{plcnt},'LineWidth',1,'MarkerSize',10)
            else
                errorbar(frac_ads(:,1,plcnt)/nch_graft,frac_ads(:,1+arch_cnt,plcnt),err_ads(:,1+arch_cnt,plcnt),...
                    'Color', pclr{arch_cnt},'Marker',msty{arch_cnt},'MarkerFaceColor',...
                    pclr{arch_cnt},'LineStyle',lsty{plcnt},'LineWidth',1,'MarkerSize',10)
            end
            legendinfo{lcnt} = [leg_arr{arch_cnt} ', PDI = ' num2str(pdi_plot(plcnt,1),'%.1f')];
            lcnt = lcnt + 1;
            
        end
        
    end
    
    legend(legendinfo,'FontSize',12,'Location','NorthWest','Interpreter','Latex')
    legend boxoff
    saveas(h1,'./../../Figs_paper/fads_npabynpc_pdi_arch.png');
    clear legendinfo
    
end

%% <Qb>
if figqnet
    
    fprintf('%s\n','Preparing net bound charge plot');
    
    % Create arrays for according to the input array sizes
    qnet_brsh = zeros(length(nfreearr),length(arch_arr)+1,length(pdi_freearr)); %+1 for y-dimension to incorporate the N_f values
    err_qnet  = zeros(length(nfreearr),length(arch_arr)+1,length(pdi_freearr)); %+1 for y-dimension to incorporate the N_f values
    pdi_plot = zeros(length(pdi_freearr));
    
    for pdi_cntr = 1:length(pdi_freearr) % begin pdi free loop
        
        pdifree     = pdi_freearr(pdi_cntr);
        pdifree_str = num2str(pdifree,'%1.1f');
        
        qnet_brsh(:,1,pdi_cntr) = nfreearr;
        err_qnet(:,1,pdi_cntr)  = nfreearr;
        pdi_plot(pdi_cntr,1)   = pdifree;
        
        % Check file existence
        fname = sprintf('./../../net_charge/AvgQnetBound_%s.txt',pdifree_str);
        fads_id = fopen(fname);
        
        if fads_id <= 0
            fprintf('%s does not exist', fname);
            continue;
        end
        
        fprintf('File under process %s\n', fname);
        
        % Read and parse first line
        tline = fgetl(fads_id); % get header
        if ~ischar(tline) || isempty(tline)
            fprintf('ERROR: Unable to read file %s\n', tline);
            continue;
        end
        spl_tline = strtrim(strsplit(strtrim(tline)));
        
        nf_col = -1; arch_col = -1; favg_col = -1; err_col = -1;ncases_col = -1;
        for wcnt = 1:length(spl_tline)
            if strcmp(spl_tline{wcnt},'N_pa')
                nf_col = wcnt;
            elseif strcmp(spl_tline{wcnt},'Arch')
                arch_col = wcnt;
            elseif strcmp(spl_tline{wcnt},'num_cases')
                ncases_col = wcnt;
            elseif strcmp(spl_tline{wcnt},'<Q_{b}>')
                favg_col = wcnt;
            elseif strcmp(spl_tline{wcnt},'StdErrMean')
                err_col = wcnt;
            end
        end
        
        if nf_col == -1 || arch_col == -1 || favg_col == -1 || err_col == -1 || ncases_col == -1
            fprintf('Could not find all headers in %s \n' , tline);
            continue
        end
        
        
        while ~feof(fads_id)
            
            % Store data into arrays by comparing with the input arrays
            tline = fgetl(fads_id); % get line
            if ~ischar(tline) && isempty(tline)
                fprintf('ERROR: Unable to read line %s\n', tline);
                continue;
            end
            spl_tline = strtrim(strsplit(strtrim(tline)));
            
            % nf value
            findnf = -1;
            for itercnt = 1:length(nfreearr) % comparing wrt the INPUT array
                if str2double(spl_tline{nf_col}) == nfreearr(itercnt)
                    findnf = 1;
                    rownum = itercnt;
                    break;
                end
            end
            
            if findnf == -1
                fprintf('WARNING: Did not find N_f value in the input array: %d\n', str2double(spl_tline{nf_col}));
                continue;
            end
            
            %arch value
            findarch = -1;
            for itercnt = 1:length(arch_arr) %comparing wrt the INPUT Array
                if strcmp(arch_arr{itercnt},spl_tline{arch_col})
                    findarch = 1;
                    colnum = itercnt;
                    break;
                end
            end
            
            if findarch == -1
                fprintf('WARNING: Did not find arch value in the input array: %s\n', spl_tline{arch_col});
                continue;
            end
            
            qnet_brsh(rownum,colnum+1,pdi_cntr) = str2double(spl_tline{favg_col});
            err_qnet(rownum,colnum+1,pdi_cntr)  = str2double(spl_tline{err_col});
            
        end
        
        fprintf('Finished analyzing %s\n', fname);
        
    end
    
    % Plot qnet, err_qnet in one go
    
    h1 = figure;
    hold on
    box on
    set(gca,'FontSize',16)
    xlabel('$n_{pa}/n_{pc}$','FontSize',20,'Interpreter','Latex')
    ylabel('$\langle Q_{b} \rangle$','FontSize',20,'Interpreter','Latex')
    
    lcnt = 1;
    for plcnt = 1:length(pdi_freearr)
        
        for arch_cnt = 1:length(arch_arr)
            
            if pdi_freearr(plcnt) == 1
                errorbar(qnet_brsh(:,1,plcnt)/nch_graft,qnet_brsh(:,1+arch_cnt,plcnt),err_qnet(:,1+arch_cnt,plcnt),...
                    'Color', pclr{arch_cnt},'Marker',msty{arch_cnt},'MarkerEdgeColor',pclr{arch_cnt},'MarkerFaceColor',...
                    'none','LineStyle',lsty{plcnt},'LineWidth',1,'MarkerSize',10)
            else
                errorbar(qnet_brsh(:,1,plcnt)/nch_graft,qnet_brsh(:,1+arch_cnt,plcnt),err_qnet(:,1+arch_cnt,plcnt),...
                    'Color', pclr{arch_cnt},'Marker',msty{arch_cnt},'MarkerFaceColor',...
                    pclr{arch_cnt},'LineStyle',lsty{plcnt},'LineWidth',1,'MarkerSize',10)
            end
            legendinfo{lcnt} = [leg_arr{arch_cnt} ', PDI = ' num2str(pdi_plot(plcnt,1),'%.1f')];
            lcnt = lcnt + 1;
            
        end
        
    end
    
    lgd = legend(legendinfo,'FontSize',12,'Location','NorthWest','Interpreter','Latex');
    legend boxoff
    saveas(h1,'./../../Figs_paper/fig_qnetbrush.png');
    
end % end if fig2b


%% \rho
if figdens
    
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
                rho_prefix = 'grpdens_config_*.lammpstrj';
                rho_fylelist = dir(strcat(dirname,'/',rho_prefix));
                if min(size(rho_fylelist)) == 0
                    fprintf('No files/Empty files are found for %s\n',rho_prefix);
                    continue;
                end
                
                nfyles = numel(rho_fylelist); %number of files of the type
                if nfyles == 0
                    fprintf('Did not find files of the type grpdens_config* in %d\t%s\n', nplot,dirstr);
                    continue;
                else
                    [latest_fyleindex] = find_latest_fyle(rho_fylelist,nfyles);
                end
                
                rho_fylename = strcat(dirname,'/',rho_fylelist(latest_fyleindex).name);
                if exist(rho_fylename,'file') ~= 2
                    fprintf('%s does not exist/empty file\n',rho_fylename);
                    continue;
                elseif struct(dir(rho_fylename)).bytes == 0
                    fprintf('Empty file: %s \n',rho_fylename);
                    continue;
                end
                fprintf('Plotting density profiles using %s for n_pa = %d\n', rho_fylename, nplot);
                
                all_data = importdata(rho_fylename,' ',1);
                fld = all_data.data;
                rdata = fld(:,1); lz = fld(1,1) + fld(length(fld(:,1)),1); %because of binning
                nbins = length(rdata(:,1));
                
                % sanity check for length of data
                if casecntr == 1
                    nbins_old = nbins;
                else
                    if nbins_old ~= nbins
                        fprintf('ERROR: The values of bins do not match %s\t%d\t%g\n',dirstr,nplot,pdifree)
                        error('CHECK nbins')
                    else
                        nbins_old = nbins;
                    end
                end
                
                avg_rho_graft = fld(:,2) + avg_rho_graft;
                avg_rho_free  = fld(:,3) + avg_rho_free;
                
                totcases = totcases + 1;
                
            end
            
            avg_rho_graft =  avg_rho_graft/totcases;
            avg_rho_free  =  avg_rho_free/totcases;
            
            plot(rdata/lz,avg_rho_graft,'Color','k','LineStyle',lsty{pdi_cntr},'LineWidth',2)
            plot(rdata/lz,avg_rho_free,'Color',green,'LineStyle',lsty{pdi_cntr},'LineWidth',2)
            
            legendinfo{lcnt} = ['Brush Monomers; ' 'PDI = ' num2str(pdifree,'%.1f')];
            lcnt = lcnt + 1;
            legendinfo{lcnt} = ['Free Monomers; ' 'PDI = ' num2str(pdifree,'%.1f')];
            lcnt = lcnt + 1;
            
        end
        
        lgd = legend(legendinfo,'FontSize',12,'Location','NorthEast','Interpreter','Latex');
        legend boxoff
        saveas(h1,sprintf('./../../Figs_paper/SI_Figs/fig_dens_%s_%d.png',dirstr,nplot));
        
    end
    
end


