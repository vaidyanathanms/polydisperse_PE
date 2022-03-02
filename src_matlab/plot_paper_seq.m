%% Plot all data for the paper

% figfads : f_ads vs Npa/Npc for different PDI and architectures
% adsorbed, entire chain is considered to be adsorbed.
% figads_mon2: If even one of 30 is adsorbed, all 30 is counted
% fignumavg_MW: number average MW based on figadsmon2 (that is even if one
% monomer is adsorbed, entire MW counts towards the calculation)

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
nfreearr = [150];
casearr  = [1,2,3,4];
pdi_freearr = [1,1.5];
pdi_strarr  = {'1.0','1.5'};
arch_arr = {'bl_bl','bl_al','al_bl','al_al'};
leg_arr  = {'Block-Block','Block-Alter','Alter-Block','Alter-Alter'}; % ALWAYS CHECK for correspondence with arch_arr
pdigraft = 1.0;
nmonfree = 30; nmongraft = 30; ngraft = 32;
cutoff = '1.50';
nch_graft = 32;
lz = 120; area=35^2;
set_tmax = 3e7; % maximum timestep for analysis;
tot_graftmon = nmongraft*ngraft;
% ONLY FOR ADSORBED CHAIN CALCULATION. TO OVERCOME THE STUPID MISTAKE IN THE ANALYSIS CODE WHERE I MISSED A FACTOR OF densfreq for adsorbed chain calculation and not for actual density calculation
densfreq = 5;

%% Input flags
% see definitions above
figfads = 0; figfads_mon2 = 0; fignumavg_MW=1;

%% Pre-calculations
rhofree = nfreearr*30/(lz*area);
pdigraft_str = num2str(pdigraft,'%1.1f');

%% f_ads (adsorbed chain fraction).
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
        fname = sprintf('./../../seq_analysis/overall/adsorbed_chain_ave_allcases_rcut_%s_pdifree_%g.dat',...
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
    xlabel('Polydispersity Index','FontSize',20,'Interpreter','Latex')
    ylabel('$f_{\rm{ads}}^{\rm{ch}}$','FontSize',20,'Interpreter','Latex')
    
    frac_combined = [frac_ads(:,2,1) frac_ads(:,3,1) frac_ads(:,4,1) frac_ads(:,5,1); 
        frac_ads(:,2,2) frac_ads(:,3,2) frac_ads(:,4,2) frac_ads(:,5,2)];
    X = categorical(pdi_strarr);
    bar(X,frac_combined)
    ylim([0 1.8]);        
    lgd = legend(leg_arr,'FontSize',16,'Location','NorthEast','Interpreter','Latex','Orientation','Horizontal');
    lgd.NumColumns = 2;

    saveas(h1,'./../../Figs_paper/fads_npabynpc_pdi_archseq.png');
    clear legendinfo
    
end
%% Adsorbed monomer with respect to the chain definition
% Definition of # of adsorbed monomers = length of the chain if one of the
% monomer of free chain is within cutoff of the graft chain. Similarly. The
% a chain is adsorbed if any of the monomer of the free is within cutoff of
% the graft chain. If 20 of 30 monomers are within cutoff, the number of adsorbed =
% 30
if figfads_mon2
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
        fname = sprintf('./../../seq_analysis/overall/adsorbedmon_wrtch_ave_allcases_rcut_%s_pdifree_%g.dat',...
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
            elseif strcmp(spl_tline{wcnt},'avg_fraction') %actually this is the number, not the fraction
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
    xlabel('Polydispersity Index','FontSize',20,'Interpreter','Latex')
    ylabel('$f_{\rm{ads}}^{\rm{mon}}$','FontSize',20,'Interpreter','Latex')
    
    frac_combined = [frac_ads(:,2,1) frac_ads(:,3,1) frac_ads(:,4,1) frac_ads(:,5,1);
        frac_ads(:,2,2) frac_ads(:,3,2) frac_ads(:,4,2) frac_ads(:,5,2)];
    X = categorical(pdi_strarr);
    bar(X,frac_combined/tot_graftmon)
    ylim([0 1.9])
    lgd = legend(leg_arr,'FontSize',16,'Location','NorthWest','Interpreter','Latex','Orientation','Horizontal');
    lgd.NumColumns = 2;
    saveas(h1,'./../../Figs_paper/fadsmon_wrtch_npabynpc_pdi_archseq.png');
    clear legendinfo
    
end

%% Plot number averaged molecular weight
% Definition of # of adsorbed monomers = length of the chain if one of the
% monomer of free chain is within cutoff of the graft chain. Similarly. The
% a chain is adsorbed if any of the monomer of the free is within cutoff of
% the graft chain. If 20 of 30 monomers are within cutoff, the number of adsorbed =
% 30
if fignumavg_MW
    
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
        fname = sprintf('./../../seq_analysis/overall/numavgMW_ave_allcases_rcut_%s_pdifree_%g.dat',...
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
            elseif strcmp(spl_tline{wcnt},'avg_fraction_numavgMW')
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
    xlabel('Polydispersity Index','FontSize',20,'Interpreter','Latex')
    ylabel('$N_{n}^{\rm ads}/N_{n}$','FontSize',20,'Interpreter','Latex')
    
    frac_combined = [frac_ads(:,2,1) frac_ads(:,3,1) frac_ads(:,4,1) frac_ads(:,5,1); 
        frac_ads(:,2,2) frac_ads(:,3,2) frac_ads(:,4,2) frac_ads(:,5,2)];
    X = categorical(pdi_strarr);
    bar(X,frac_combined/nmonfree)
    ylim([0 2]);        
    lgd = legend(leg_arr,'FontSize',16,'Location','NorthWest','Interpreter','Latex');
    
    saveas(h1,'./../../Figs_paper/fnumavgMW_npabynpc_pdi_seq.png');
end

