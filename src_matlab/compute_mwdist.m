function [distout] = compute_mwdist(molarr,column_num)

% To compute the molecular weight distribution.
% Inputs: array name and the column number for which distn needs to be
% calculated.

% Keep bin width to 8.
% Filter out ZERO molecular weight

[~,~,filt_dist] = find(molarr(:,column_num));
binwid    = 8;
distout   = histogram(filt_dist,'BinWidth',binwid,'Normalization','pdf');