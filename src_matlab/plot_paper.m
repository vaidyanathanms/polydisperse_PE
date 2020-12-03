%% Plot all data for the paper

% Fig 1  : SZ-distribution
% Fig 2a : f_ads vs Npa/Npc for different PDI and architectures
% Fig 2b : Qnet vs Npa/Npc for different PDI and architectures
% Fig 3a : p_ads(MW) vs MW for Block-Block and various N_pa/Npc
% Fig 3b : p_ads(MW) vs MW for Alter-Alter and various N_pa/Npc
% Fig 4  : Density plots
% Fig 5  : Bidispersed case

clear
clc
close all
format long

%% Color codes
green = [0 0.5 0.0]; gold = [0.9 0.75 0]; orange = [0.91 0.41 0.17];brown = [0.2 0 0];
pclr = {'r',green,'b',brown,'k','m', gold};
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
fig1  = 0;
fig2a = 1;
fig2b = 1;
fig3a = 0;
fig3b = 0;
fig4  = 0;

%% Pre-calculations
rhofree = nfreearr*30/(lz*area);
pdigraft_str = num2str(pdigraft,'%1.1f');

%% Fig1
if fig1
    
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

%% Fig2a
if fig2a
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

%% Fig2b
if fig2b
    
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


%% Fig3a
if fig3a
    %Run avg_and_coarsen_dist.m in conjunction with out_mwdist.m
end

%% Fig3b
if fig3b
    %Run avg_and_coarsen_dist.m in conjunction with out_mwdist.m
end

%% Fig 4
if fig4
   % Run density_plot.m 
end
   
%% Fig 5
if fig4
    % Run bidispersed_analysis.m
end

