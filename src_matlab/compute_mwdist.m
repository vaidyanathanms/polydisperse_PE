function [distout] = compute_mwdist(molarr,column_num)

% To compute the molecular weight distribution.
% Inputs: array name and the column number for which distn needs to be
% calculated.

% Arbitrarily keep num of bins at one fourth the max mol. wt. 
% So for uniform dist, every 4 mol.wt will contribute to one bin.
% Filter out ZERO molecular weight

[~,~,filt_dist] = find(molarr(:,column_num))
numbins   = floor(0.25*max(filt_dist))
distout   = histogram(filt_dist,numbins);



