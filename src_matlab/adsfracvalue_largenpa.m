%% To create histograms for chain MW

clc;
clear;
close all;
format long

%% Input data

nfree = 150;
orange = [0.91 0.41 0.17]; brown=[0.6 0.2 0];  gold = [0.9 0.75 0];

pdixvals = [1 1.5];
fracvals = [1.4578 1.3314; 0.82719	1.0778];


hpdi = figure;
hold on
box on
set(gca,'FontSize',16)
xlabel('Polydispersity Index','FontSize',16,'Interpreter','Latex')
ylabel('$f_{ads}$','FontSize',16,'Interpreter','Latex')
bar(pdixvals,fracvals)
legend('Block/Block','Alter/Alter')
legend boxoff
saveas(hpdi,sprintf('fads_n_%d',nfree),'png');