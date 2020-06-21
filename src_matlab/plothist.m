%% To plot various histograms of output data

clc;
close all;
clear;
format long

%% Color codes
green = [0 0.5 0.0]; gold = [0.9 0.75 0]; orange = [0.91 0.41 0.17];brown = [0.2 0 0];
pclr = {'m',brown,green,'k','m', gold};
lsty = {'-','--',':'};
msty = {'d','s','o','x'};

%% Inputs
maxnum_cases = 4;
pdi_freearr = [1,1.3,1.5];
pdigraft = 1.0;
nmonfree = 30; nmongraft = 30; ngraft = 32;
cutoff = '1.50';
lz = 120; area=35^2;
set_tmax = 3e7; % maximum timestep for analysis;

%% Input flags

comppdi   = 1; % compare 2 pdis for same arch.
comparch  = 1; % compare 2 archs for same pdi.

% Arch names DO NOT have an underscore separating the graft and free arch.

pdi_comp_cell  = {'1.0','1.5','blbl'}; % 2 pdis for comparison followed by arch
arch_comp_cell = {'blbl','alal','1.0'}; % 2 archs for comparison followed by pdi

%% Pre-calculations

pdigraft_str = num2str(pdigraft,'%1.1f');

%% Main Analysis

if comppdi % compare pdi fixing architecture
    
    fylename  = sprintf('./../../All_analysis.xlsx');
    if ~isfile(fylename)
        error('%s does not exist', fylename);
    end
    sheetname = sprintf('compare_pdi_%s_%s_%s',pdi_comp_cell{1},pdi_comp_cell{2},pdi_comp_cell{3});
    maindata  = readtable(fylename,'Sheet',sheetname,'Range','A:I','ReadVariableNames',1);
    
end
