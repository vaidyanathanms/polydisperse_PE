% Plot SZ distribution

clc;
clear;
close all;
format long

mnarr = 10:10:30;
pdifree = 2;
figure
hold on
box on
for i = 1:length(mnarr)
    mnideal = mnarr(i);
    
    mwvals = 1:10*mnideal;
    mwrat  = mwvals/mnideal;
    k = 1/(pdifree - 1);
    
    term1 = k^k;
    term2 = mwrat.^(k-1);
    term3 = exp(-k*mwrat);
    term4 = gamma(k);
    
    psztheory = (term1.*term2.*term3)./term4;
    normval = trapz(mwvals,psztheory);
    
    plot(mwvals,psztheory/normval,'LineWidth',2);
    legendinfo{i} = num2str(mnideal);
end
xlim([0 800])
legend(legendinfo,'FontSize',16,'Location','Best','Interpreter','Latex')
clear legendinfo

pdiarr = 1.5:1:1.5;
mnideal = 4000;
mwvals = 1:10*mnideal;
mwrat  = mwvals/mnideal;
figure
hold on
box on
for j = 1:length(pdiarr)
    pdifree = pdiarr(j);
    k = 1/(pdifree - 1);
    
    term1 = k^k;
    term2 = mwrat.^(k-1);
    term3 = exp(-k*mwrat);
    term4 = gamma(k);
    
    psztheory2 = (term1.*term2.*term3)./term4;
    normval = trapz(mwvals,psztheory2);
    
    plot(mwvals,psztheory2/normval,'LineWidth',2);
    legendinfo{j} = ['PDI: ' num2str(pdifree)];
end
set(gca,'xscale','log')
xlim([10 40000])
set(gca,'FontSize',16)
xlabel('Mol. Wt','FontSize',20,'Interpreter','Latex')
ylabel('Probability','FontSize',20,'Interpreter','Latex')
legend(legendinfo,'FontSize',16,'Location','Best','Interpreter','Latex')
