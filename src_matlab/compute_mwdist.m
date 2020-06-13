function [distout] = compute_mwdist(molarr)

% To compute the molecular weight distribution.

% Arbitrarily keep num of bins at one fourth the max mol. wt. 
% So for uniform dist, every 4 mol.wt will contribute to one bin.
maxbinval = floor(0.25*molarr(:,3)); 
distout = histogram(molarr(:,3),maxbinval);



