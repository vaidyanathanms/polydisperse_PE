function [molcnt] = analyze_datafile(dataname,nfree_ch)

% Check whether file exists
fid = fopen(dataname,'r');

if fid <= 0
    fprintf('%s not found',dataname);
    return
end

fprintf('Analyzing %s\n', dataname);
natom_flag = -1; % to see natoms keyword is read before reading Atoms
molcnt = zeros(nfree_ch,1);

% Start reading file
fgetl(fid); % skip first line;

while true
    
    tline = fgetl(fid);
    if ~ischar(tline) || isempty(tline)
        continue;
    end
    spl_tline = strsplit(strtrim(tline));
        
    if length(spl_tline) > 1 % for double keywords check this first
        if strcmp(spl_tline{2},'atoms') % check for atoms keyword
        
            num_atoms = str2double(spl_tline{1});
            natom_flag = 1;
            
        end
        
    elseif strcmp(spl_tline{1},'Atoms') % check for Atoms keyword
        
        if natom_flag == -1
            fprintf('%s Did not find natoms keyword in', dataname);
            return
        end
        
        for i = 1:num_atoms % increment molcnt counter based on molid
            tline = fgetl(fid); spl_tline = strsplit(strtrim(tline));
            molid = str2double(spl_tline{2});
            molcnt(molid) = molcnt(molid) + 1;
        end
        
        break; % no need to read the rest of the file
        
    end
    
end
