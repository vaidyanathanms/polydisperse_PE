function [norm_distcounts,sorted_init_uniq_arr] = find_distribution_of_mw(dataarr,init_dist,total_frames)

%% Format
% Input
% dataarr -> 1D array containing the MW of all chains that are
% adsorbed over the course of the simulation (summing across cases and
% num_fyles per case for a given arch).
% init_dist: 1D array acontaining the MW of all chains at the beginning of
% the simulation (summing across all the cases for a given arch).
% total_frames: total number of frames analyzed (summing across all the
% cases and num_fyles per case for a given arch.)

%Outputs: norm_distcounts (4D array): col1: list of all initial MWs.
%identical to the first column of sorted_init_uniq_arr; col2: total
%occurence of a given mol wt; col3: total occurence per frame; col4: total
%occurence/frame/inital number of chains with a given MW
%sorted_init_uniq_arr: (2D array): col1: list of all initial MWs, col2:
%occurence of each MW at the beginning. sum(col2)=num chains
%% Analysis

% Counts of initial MW distribution: https://www.mathworks.com/help/matlab/ref/unique.html
[init_unique,~,init_ic] = unique(init_dist(:,1));
all_init_counts = accumarray(init_ic,1);
init_unique_counts = [init_unique, all_init_counts];

% Counts of adsorbed MW distribution
dataarr(:,1) = dataarr(dataarr(:,1)~=0,1); % weed out zero values
[adsorbed_unique,~,ads_ic] = unique(dataarr(:,1));
all_adsorbed_counts = accumarray(ads_ic,1);
adsorbed_unique_counts = [adsorbed_unique, all_adsorbed_counts];

% Sort based on the MW from smallest to largest for both init and adsorbed.
sorted_init_uniq_arr = sortrows(init_unique_counts,1);
sorted_adsb_uniq_arr = sortrows(adsorbed_unique_counts,1);

% Run through all the elements in sorted_adsb_uniq_arr and divide by the
% corresponding counts in the sorted_init_uniq_arr to compute the
% occurences per MW. Then divide by the number of frames to normalize.
% p(n_ads^mw) = [\sum_{i=1}^{n_f^mw} n_i,ads^MW \times num_frames]/[n_f^MW
% \times \sum N_f]
norm_distcounts = zeros(length(sorted_init_uniq_arr),3);
norm_distcounts(:,1) = sorted_init_uniq_arr(:,1);

for i = 1:length(sorted_adsb_uniq_arr)
    adsorb_MW = sorted_adsb_uniq_arr(i,1);
    adsorbed_num_chains = sorted_adsb_uniq_arr(i,2);
    find_in_init = -1;
    
    for j = 1:length(sorted_init_uniq_arr) %find which element in init array does this correspond to
        init_MW = sorted_init_uniq_arr(j,1);
        
        if adsorb_MW == init_MW
            find_in_init = 1;
            init_index = j;
            init_num_chains = sorted_init_uniq_arr(j,2); % number of chains at the beginning with the given adsorbed MW
            break; % no need to search the rest
        end
        
    end
    
    if find_in_init == -1
        fprintf('ERROR: the molecular wt of adsorbed chain (%d) not found in the initial MW array',adsorb_MW);
        return;
    end
    
    norm_distcounts(init_index,2) = adsorbed_num_chains; %total occurence
    norm_distcounts(init_index,3) = adsorbed_num_chains/total_frames; %total occurence per frame
    norm_distcounts(init_index,4) = adsorbed_num_chains/(init_num_chains*total_frames); %total occurence/frame/initial number of same MW
    
end

