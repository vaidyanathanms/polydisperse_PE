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
pclr = {'m',brown,green,'k','m', gold};
lsty = {'-','--',':'};
msty = {'d','s','o','x'};

%% Input data
nfreearr = [16;32;64;96;128;150];
casearr  = [1,2,3,4];
pdi_freearr = [1,1.5];
arch_arr = {'bl_bl','al_al'};
leg_arr  = {'Block-Block','Alter-Alter'}; % ALWAYS CHECK for correspondence with arch_arr
pdigraft = 1.0;
nmonfree = 30; nmongraft = 30; ngraft = 32;
cutoff = '1.50';
lz = 120; area=35^2;
set_tmax = 3e7; % maximum timestep for analysis;

%% Input flags
% see definitions above
fig1  = 1; 
fig2a = 0; 
fig2b = 0; 
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
    
end

%% Fig2b
if fig2b
    
end

%% Fig3a
if fig3a
    
end

%% Fig3b
if fig3b
    
end

%% Fig 4
if fig4
    
end