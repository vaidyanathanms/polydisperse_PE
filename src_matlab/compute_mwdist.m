function [distout] = compute_mwdist(molarr,column_num)

% To compute the molecular weight distribution.
% Inputs: array name and the column number for which distn needs to be
% calculated.

% Arbitrarily keep num of bins at one fourth the max mol. wt. 
% So for uniform dist, every 4 mol.wt will contribute to one bin.
maxbinval = floor(0.25*max(molarr(:,column_num))); 
distout   = histogram(molarr(:,column_num),maxbinval);



