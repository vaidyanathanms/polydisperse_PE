function nmonlist = find_all_mons_name(dataname,idinplist)
% Required order in datafile: graft mons -> free mons -> counter ions -> salt

nmonlist = zeros(length(idinplist),1);
% Check whether file exists
fid = fopen(dataname,'r');

if fid <= 0
    fprintf('%s not found',dataname);
    return
end

% Initialize data
fprintf('Analyzing datafile: %s\n', dataname);
natom_flag = -1; % to see natoms keyword is read before reading Atoms

% Start reading file
fgetl(fid); % skip first line;

while true
    
    tline = fgetl(fid);
    if ~ischar(tline) || isempty(tline)
        continue;
    end
    spl_tline = strsplit(strtrim(tline));
    
    if strcmp(spl_tline{1},'Atoms') % check for Atoms keyword
        
        if natom_flag == -1
            fprintf('%s Did not find natoms keyword in', dataname);
            return
        end
        
        i = 1;
        while (i <= num_atoms) % increment molcnt counter based on molid
            
            if feof(fid)
                break;
            end
            
            tline = fgetl(fid);
            spl_tline = strsplit(strtrim(tline));
            if isempty(tline) || ~ischar(tline)
                continue;
            end
            
            atype = str2double(spl_tline{3});
            for u = 1:length(idinplist)
                if atype == idinplist(u)
                    nmonlist(u,1) = nmonlist(u,1) + 1;
                    break;
                end
            end
            i = i + 1;
            
        end
        
        break; % no need to read the rest of the file
        
    elseif length(spl_tline) > 1 % for double keywords check this first
        
        if strcmp(spl_tline{2},'atoms') % check for atoms keyword
            
            num_atoms = str2double(spl_tline{1});
            natom_flag = 1;
            
        end % do this "ELSEIF" condition after the first "IF" condition, so that the elseif statement works correctly.
        
    end
    
end
