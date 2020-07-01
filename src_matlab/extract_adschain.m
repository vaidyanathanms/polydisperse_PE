function [adschain_arr] = extract_adschain(main_fylename,max_mw_ads)

% Required header in main_fyle: "freechainMW" with comma separated 
% ONLY HEADER IS COMMA DELIMITED

% Check whether file exists
fid = fopen(main_fylename,'r');

if fid <= 0
    fprintf('%s not found \n',main_fylename);
    return
end

% Initialize data
fprintf('Analyzing %s\n', main_fylename);
adschain_arr = zeros(max_mw_ads,2); % Maximum value of the MW that can be adsorbed
adschain_arr(:,1) = 1:max_mw_ads; % Transpose is automatically taken care of by MATLAB

% Start reading file and find header keyword
tline = fgetl(fid); % get header
if ~ischar(tline) || isempty(tline)
    fprintf('ERROR: Unable to read file %s\n', tline)
    return;
end
spl_tline = strsplit(strtrim(tline),','); % comma demilited
find_keyword = -1; % to find "freechainMW" keyword

for wordcnt = 1:length(spl_tline)
    if strcmp(spl_tline{wordcnt},'freechainMW')
        find_keyword = 1; column_num = wordcnt;
        break;
    end
end
if find_keyword == -1
    fprintf('ERROR: Unable to find keyword (freechainMW) in header %s\n', tline);
    return;
end
    
% start reading the rest of the file
while true
    
    tline = fgetl(fid); % get header
    if ~ischar(tline) || isempty(tline)
        continue;
    end
    
    spl_tline = strsplit(strtrim(tline)); 
    ch_mw = str2double(spl_tline{column_num});
    
    % check if this exceeds the max designated value
    if ch_mw > max_mw_ads 
        fprintf('ERROR: MW of ads_chains (%d) exceeds the max value (%d). Change "max_mw_free" value in the "Pre-calculations" section',...
            ch_mw, max_mw_ads);
        break;
    end
 
    % check if zero/negative or real values are there
    if ch_mw < 1 || ~isinteger(ch_mw)
        fprintf('ERROR: Unknown output MW %d\n',ch_mw)
        break;
    end
    
    adschain_arr(ch_mw,2) = adschain_arr(ch_mw,2) + 1; % add corresponding value
    
end
    
    
    
    