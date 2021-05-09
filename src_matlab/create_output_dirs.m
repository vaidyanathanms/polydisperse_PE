function [stat] = create_output_dirs(dirname)
% Create output directories
stat = 1; % directory exists
if ~exist(dirname,'dir')
    fprintf('Making %s\n', dirname);
    stat = mkdir(dirname);
end