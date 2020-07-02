function [distout] = compute_mwdist(dataarr,column_num)

% To compute the molecular weight distribution.
% Inputs: array name and the column number for which distn needs to be
% calculated.

% Keep bin width to 8.
% * BUG FIX: DO NOT EDIT OUT ZERO MOLECULAR WT
% Earlier version had zero MW wt filtered out. That doesn't make sense. If
% there are no MWs between 0 and 8, the first bin has to be identically
% zero. Good way os to add the bin limits. 

binwid    = 8;
distout   = histogram(dataarr(:,column_num),'BinWidth',binwid,'Normalization','pdf','BinLimits',[1,max(filt_dist)+1]);