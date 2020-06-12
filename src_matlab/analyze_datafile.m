function [molcnt] = analyze_datafile(dataname,nfree_ch)

% Check whether file exists
fid = fopen(dataname,'r');

if fid <= 0
    fprintf('%s not found',dataname);
    return
end

natom_flag = -1; % to see natoms keyword is read before reading Atoms
molcnt = zeros(nfree_ch,1); 

% Start reading file
while fid < 0
    tline = fgetl(finp);
    if ~ischar(tline)
        continue;
    end
    spl_tline = strsplit(strtrim(tline));
    
    if strcmp(spl_tline{2},'natoms') % check for natoms keyword
        
        num_atoms = str2double(spl_tline{1});
        natom_flag = 1;
        
    elseif strcmp(spl_tline{1},'Atoms') % check for Atoms keyword
        
        if natom_flag == -1
            fprintf('%s Did not find natoms keyword in', dataname);
            return
        end
        
        for i = 1:num_atoms % increment molcnt counter based on molid
            molid = str2double(spl_tline{2});
            molcnt(molid) = molcnt(molid) + 1;
        end 
        
    end
        
end
