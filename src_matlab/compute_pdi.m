function [pdival,mn_val,mw_val] = compute_pdi(molarr,nfree_ch)

% https://en.wikipedia.org/wiki/Dispersity
% PDI = (\Sigma Ni * \Sigma Ni Mi^2)/(\Sigma Ni Mi)^2

% check whether any of the third column has zero value - sanity check

if ~isempty(find(~molarr,3))
    fprintf('Found zero molecular weight values in molarr \n');
    return;
end

sig_ni = nfree_ch;
sig_nimi = sum(molarr(:,3));
sig_nimisq = sum((molarr(:,3).^2));

mn_val = sig_nimi/sig_ni;
mw_val = sig_nimisq/sig_nimi;
pdival = mw_val/mn_val;
    