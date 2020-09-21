%% Plot all data for the paper

% Fig 1: SZ-distribution
% Fig 2a : f_ads vs Npa/Npc for different PDI and architectures
% Fig 2b : Qnet vs Npa/Npc for different PDI and architectures
% Fig 3a : p_ads(MW) vs MW for Block-Block and various N_pa/Npc
% Fig 3b : p_ads(MW) vs MW for Alter-Alter and various N_pa/Npc
% Fig 4: Bidispersed case

clear
clc
close all
format long

%% Color codes
green = [0 0.5 0.0]; gold = [0.9 0.75 0]; orange = [0.91 0.41 0.17];brown = [0.2 0 0];
pclr = {'r',green,'b',brown,'k','m', gold};
lsty = {'-','--',':'};
msty = {'d','o','s','x'};

%% Input data
nfreearr = [16]%;32;64;128;150];
casearr  = [1]%,2,3,4];
pdi_freearr = [1,1.5];
arch_arr = {'bl_bl'}%,'al_al'};
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
fig2a = 0;
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
        legendinfo{plcnt} = ['$N_{pa} =$ ' num2str(nval_pl(plcnt))];
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
    xlabel('$N_{pa}/N_{pc}$','FontSize',20,'Interpreter','Latex')
    ylabel('$f_{\rm{ads}}$','FontSize',20,'Interpreter','Latex')
    
    lcnt = 1;
    for plcnt = 1:length(pdi_freearr)
        
        for arch_cnt = 1:length(arch_arr)
            
            errorbar(frac_ads(:,1,plcnt)/nch_graft,frac_ads(:,1+arch_cnt,plcnt),err_ads(:,1+arch_cnt,plcnt),...
                'Color', pclr{arch_cnt},'Marker',msty{arch_cnt},'MarkerFaceColor',...
                pclr{arch_cnt},'LineStyle',lsty{plcnt},'LineWidth',2,'MarkerSize',8)
            
            legendinfo{lcnt} = [leg_arr{arch_cnt} ', PDI = ' num2str(pdi_plot(plcnt,1))];
            lcnt = lcnt + 1;
            
        end
        
    end
    
    legend(legendinfo,'FontSize',12,'Location','NorthWest','Interpreter','Latex')
    legend boxoff
    saveas(h1,'./../../Figs_paper/fads_npabynpc_pdi_arch.png');
    
end

%% Fig2b 
if fig2b
    
    fprintf('%s\n','Preparing net bound charge plot');
    chargearr = [0;0;1;0;-1;1;-1];
    rbin = 1;

    h1 = figure;
    hold on
    box on
    set(gca,'FontSize',16)
    xlabel('$r$','FontSize',20,'Interpreter','Latex')
    ylabel('$f$','FontSize',20,'Interpreter','Latex')
    
    
    for pdi_cntr = 1:length(pdi_freearr) % begin pdi free loop
        
        ref_pdifree     = pdi_freearr(pdi_cntr);
        pdifree_str     = num2str(ref_pdifree,'%1.1f');
        fnr   = fopen(sprintf('./../../net_charge/QnetBound_%s.txt',pdifree_str),'w');
        fmon  = fopen(sprintf('./../../net_charge/nmonlist_%s.txt',pdifree_str),'w');
        fave  = fopen(sprintf('./../../net_charge/AvgQnetBound_%s.txt',pdifree_str),'w');
        fprintf(fnr,'%s \n','NetCharge: Q_{b}=\Delta(n_g)*\int(\sum(q_j n_j(z)dz, j=all entities)z=0,Lz))');
        fprintf(fnr,'%s\t %s\t %s\t %s\t %s\n','N_pa','Arch','Case','First tstep','Q_{b}');
        fprintf(fmon,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n',...
            'N_pa','Arch','Case','N_bb','N_ntbr','N_chbr','N_ntfr','N_chfr','Nposions','Nnegions');      
        fprintf(fave,'%s\t %s\t %s\t %s\n','N_pa','Arch','<Q_{b}>','StdDev');      
        pdigraft_str = num2str(pdigraft,'%1.1f');

        for ncnt = 1:length(nfreearr) % begin nfree loop
            nval = nfreearr(ncnt);
            
            for arch_cnt = 1:length(arch_arr)  % begin arch loop
                dirstr = arch_arr{arch_cnt};
                
                for casecntr = 1:length(casearr) % begin case loop
                    casenum = casearr(casecntr);
                
                    fprintf('Analyzing pdi/n_pa/arch/case#: %g\t%d\t%s\t%d\n',ref_pdifree,nval,dirstr,casenum)
                    
                    % Check input file name to count number of monomers of
                    % each type
                    fprintf('Extracting monomer details ..\n')
                    inp_prefix = sprintf('PEinitdata.txt');
                    inp_fylelist = dir(strcat(dirname,'/',inp_prefix));
                    if min(size(inp_fylelist)) == 0
                        fprintf('No files/Empty files are found for %s\n',inp_prefix);
                        continue;
                    end                    
                    nfyles = numel(inp_fylelist); %number of files of the type
                    if nfyles ~= 1
                        fprintf('WARNING: Found %d initial datafiles', nfyles);
                        continue;
                    end
                    % begin running through the datafile
                    inp_fylename = strcat(dirname,'/',inp_fylelist(1).name);
                    if exist(inp_fylename,'file') ~= 2
                        fprintf('%s does not exist/empty file\n',inp_fylename);
                        continue;
                    elseif struct(dir(inp_fylename)).bytes == 0
                        fprintf('Empty file: %s \n',inp_fylename);
                        continue;
                    end
                    idlist = [2;3;4;5;6;7;8]; % [nbb_brush,nneut_brush,ncharg_brush,nneut_free,ncharg_free,nposions,nnegions]...
                    outmonlist = find_all_mons_name(inp_fylename,idlist);
                    fprintf(fmon,'%d\t%s\t%d\t',nval,dirstr,casenum);
                    for u = 1:length(outmonlist)
                        fprintf(fmon,'%d\t',outmonlist(u));
                    end
                    fprintf(fmon,'\n');
                    
                    % Now analyze all the charge details
                    dirname = sprintf('./../../sim_results/outresults_dir_n_%d/%s/pdifree_%s_pdigraft_%s/Case_%d',...
                        nval,dirstr,pdifree_str,pdigraft_str,casenum);                   
                    if ~exist(dirname,'dir')
                        fprintf('%s does not exist\n',dirname);
                        continue
                    end
                    fprintf('Analyzing pdi/nfree/arch/case: %g\t%d\t%s\t%d\n',ref_pdifree,nval,dirstr,casenum);
                    
                    %check if file type exists
                    qnet_prefix = 'dens_config_*.lammpstrj';
                    qnet_fylelist = dir(strcat(dirname,'/',qnet_prefix));
                    if min(size(qnet_fylelist)) == 0
                        fprintf('No files/Empty files are found for %s\n',qnet_prefix);
                        continue;
                    end
                    
                    nfyles = numel(qnet_fylelist); %number of files of the type
                    if nfyles == 0
                        fprintf('Did not find files of the type dens_config* in %d\t%s\n', nfree,dirstr);
                    else
                        [latest_fyleindex] = find_latest_fyle(qnet_fylelist,nfyles);
                    end

                    qnet_fylename = strcat(dirname,'/',qnet_fylelist(latest_fyleindex).name);
                    if exist(qnet_fylename,'file') ~= 2
                        fprintf('%s does not exist/empty file\n',qnet_fylename);
                        continue;
                    elseif struct(dir(qnet_fylename)).bytes == 0
                        fprintf('Empty file: %s \n',qnet_fylename);
                        continue;
                    end
                    fprintf('Analyzing %s\n', qnet_fylename);
                    
                    
                    all_data = importdata(qnet_fylename);
                    fld = all_data.data;
                    rdata = fld(:,1); lz = fld(:,1) + fld(len(fld(:,1)),1); %because of binning
                    
                    ntotarr = nch_graft + sum(outmonlist);
                    pbackbon_brush = fld(:,2)*outmonlist(:,1)*lz*area/nbins;
                    pneutral_brush = fld(:,3)*outmonlist(:,2)*lz*area/nbins;
                    pnegativ_brush = fld(:,4)*outmonlist(:,3)*lz*area/nbins;
                    pneutral_free  = fld(:,5)*outmonlist(:,4)*lz*area/nbins;
                    ppositive_free = fld(:,6)*outmonlist(:,5)*lz*area/nbins;
                    pposions  = fld(:,7)*outmonlist(:,6)*lz*area/nbins;
                    pnegions  = fld(:,8)*outmonlist(:,7)*lz*area/nbins;

                    netcharge = chargearr(1)*pbackbon_brush + chargearr(2)*pneutral_brush ...
                        + chargearr(3)*pnegativ_brush + chargearr(4)*pneutral_free + ...
                        chargearr(5)*ppositive_free + chargearr(6)*pposions + chargearr(7)*pnegions;
                    
                    sumq = 0.5*netcharge(1);
                    
                    fid_g  = fopen(sprintf('./../../results_dens/results_%d_%s/grpdens.lammpstrj',nfree,dirstr));
                    data_g = textscan(fid_g,'%f%f%f','Headerlines',1);
                    
                    fld_g   = cell2mat(data_g);
                    dens_g  = fld_g(:,2);
                    [maxden, imaxden]  = max(dens_g);
                    
                    % Find edge of brush
                    for k = imaxden:length(dens_g)
                        if dens_g(k,1) < 0.05*maxden
                            i_edge = k;
                            break;
                        end
                    end
                    
                    qofr = zeros(i_edge,1); rqofr = zeros(i_edge,1);
                    qofr(1,1) = sumq; rqofr(1,1) = 0.5*rdata(1);
                    
                    % Integrate charge to the edge of the brush
                    for k = 2:i_edge
                        
                        sumq = sumq + 0.5*(rdata(k)-rdata(k-1))*(netcharge(k)+netcharge(k-1));
                        qofr(k,1)  = sumq;
                        rqofr(k,1) = 0.5*(rdata(k)+rdata(k+1));
                        
                    end
                    
                    fprintf(fnr,'%s\t%g\n',dirstr,sumq);
                    %fprintf('%s\t%g\n',dirstr,sumq);
                    fclose(fid_g);
                    fclose(fid);
                    
                    if(nvals == 1)
                        
                        plot(rqofr,qofr,'color',pclr{i},'LineWidth',2,'LineStyle',lsty{1})
                        
                    end
                    
                    
                    
                end % end case loop
                
            end % end arch loop
            
        end % end ncnt loop
        
    end % end pdi loop
    
end % end if fig2b



%                     
%                     

%                     
%                     
%             
%         end
%         
%         fclose(fnr);
%         
%     end
%     
%     % Plot Q_{b} data
%     
%     data_bb = zeros(length(nfreearr),1);
%     data_ab = zeros(length(nfreearr),1);
%     data_ba = zeros(length(nfreearr),1);
%     data_aa = zeros(length(nfreearr),1);
%     
%     for i = 1:length(nfreearr)
%         
%         fnr = fopen(sprintf('./../../all_txtfiles/fig6b_QnetBound_%d.txt',nfreearr(i)),'r');
%         data = textscan(fnr,'%s%f','Headerlines',2);
%         
%         fld = cell2mat(data(2));
%         data_bb(i,1) = fld(1,1);
%         data_ba(i,1) = fld(2,1);
%         data_ab(i,1) = fld(3,1);
%         data_aa(i,1) = fld(4,1);
%         
%         fclose(fnr);
%         
%     end
%     
%     h1 = figure;
%     hold on
%     box on
%     set(gca,'FontSize',16)
%     xlabel('$N_{pa}/N_{pc}$','FontSize',20,'Interpreter','Latex')
%     ylabel('$Q_{b}$','FontSize',20,'Interpreter','Latex')
%     
%     plot(nfreearr/ngraft,data_bb,'color',pclr{1},'LineWidth',2,'LineStyle',lsty{3},'Marker',msty{1},'MarkerSize',8,'MarkerFaceColor',pclr{1})
%     plot(nfreearr/ngraft,data_ba,'color',pclr{2},'LineWidth',2,'LineStyle',lsty{3},'Marker',msty{2},'MarkerSize',8,'MarkerFaceColor',pclr{2})
%     plot(nfreearr/ngraft,data_ab,'color',pclr{3},'LineWidth',2,'LineStyle',lsty{3},'Marker',msty{3},'MarkerSize',8,'MarkerFaceColor',pclr{3})
%     plot(nfreearr/ngraft,data_aa,'color',pclr{4},'LineWidth',2,'LineStyle',lsty{3},'Marker',msty{4},'MarkerSize',8,'MarkerFaceColor',pclr{4})
%     
%     
%     legendinfo{1} = 'Block-Block';
%     legendinfo{2} = 'Block-Alter';
%     legendinfo{3} = 'Alter-Block';
%     legendinfo{4} = 'Alter-Alter';
%     
%     legend(legendinfo,'Interpreter','Latex','FontSize',16,'Location','Best')
%     legend boxoff
%     saveas(h1,'./../../all_figures/Fig1b_QnetBound_%s','png');
%     clear legendinfo
% end

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
    
end