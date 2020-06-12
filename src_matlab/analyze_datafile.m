function [molcnt_arr] = analyze_datafile(dataname,nfree_ch,ngraft_ch)

% Order in datafile: graft mons -> free mons -> counter ions -> salt

% Check whether file exists
fid = fopen(dataname,'r');

if fid <= 0
    fprintf('%s not found',dataname);
    return
end

% Initialize data
fprintf('Analyzing %s\n', dataname);
natom_flag = -1; % to see natoms keyword is read before reading Atoms
molcnt_arr = zeros(nfree_ch,3);
max_chain_id  = ngraft_ch + nfree_ch; % maximum range of mol list from file read

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
            
            tline = fgetl(fid);
            spl_tline = strsplit(strtrim(tline));
            if ~ischar(tline) || isempty(tline)
                continue;
            end
            
            if fix(str2double(spl_tline{2})) > ngraft_ch && ...
                    fix(str2double(spl_tline{2})) <= max_chain_id %look at only free mons
                
                molid = str2double(spl_tline{2});
                remap_molid = molid - ngraft_ch;
                
                molcnt_arr(remap_molid,1) = molid; %actual molid in datafile
                molcnt_arr(remap_molid,2) = remap_molid; %remapped id
                molcnt_arr(remap_molid,3) = molcnt_arr(remap_molid,3) + 1; %mw for a particular chain
                
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
