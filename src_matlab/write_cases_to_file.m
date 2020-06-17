function outflag = write_cases_to_file(fylename,nref,pdiref,arch_ref)

% the consolidated file will be read and broken down into smaller files
% with the 
fin_cons = fopen(fylename,'r');

outflag = -1;
if fin_cons <= 0 % check for average list
    fprintf('%s does not exist', fylename);
    return;
end

headerflag = -1; % find headerflag
pdi_col = -1; nf_col = -1; arch_col = -1; casenum_col = -1;
numsam_col = -1; avgf_col = -1;

outflag = 1;
while true
    
    if headerflag == -1
        
        tline = fgetl(fout_cons);
        if ~ischar(tline) || isempty(tline)
            continue;
        end
        
        spl_tline = strsplit(strtrim(tline));
        len_tline = length(spl_tline);
        
        if len_tline ~= 8
            fprintf('Unknown number of keywords in %s \n', tline);
            continue;
        end
        
        for word_cnt = 1:len_tline % add more if keywords are added
            
            if strcmp(spl_tline{word_cnt},'PDI_free') % for pdi_free keyword
                pdi_col = word_cnt;
            elseif strcmp(spl_tline{word_cnt},'N_f') % for N_f keyword
                nf_col = word_cnt;
            elseif strcmp(spl_tline{word_cnt},'Arch') % for Arch keyword
                arch_col = word_cnt;
            elseif strcmp(spl_tline{word_cnt},'Case_num') % for Case_num keyword
                casenum_col = word_cnt;
            elseif strcmp(spl_tline{word_cnt},'numsample_pts') % for number of samples keyword
                numsam_col = word_cnt;
            elseif strcmp(spl_tline{word_cnt},'avg_fraction') % for avg_fraction keyword
                avgf_col = word_cnt;
            end
            
        end
        
        headerflag = 1; % found headerflag
        
    else
        
        tline = fgetl(fout_cons);
        spl_tline = strsplit(strtrim(tline));
        if ~ischar(tline) || isempty(tline)
            continue;
        end
        
        spl_tline = strsplit(strtrim(tline));
        len_tline = length(spl_tline);
        
        if len_tline ~= 8 && len_tline ~= 9
            fprintf('Unknown number of keywords in %s \n', tline);
            continue;
        end
        
        % parse elements
        nval   = str2double(spl_tline{nf_col});
        dirstr = spl_tline{arch_col};
        pdival = str2double(spl_tline{pdi_col});
        
        % cross check with the main array
        if ismember(nval, nfreearr)
            nval_index = find(nval==nfreearr);
        else
            fprintf('did not find %d in nfreearr\n', nval_index)
            continue;
        end
        
        if contains(arch_arr,dirstr)
            arch_index = find(contains(arch_arr,dirstr));
        else
            fprintf('did not find %s in arch_arr\n', dirstr)
            continue;
        end
        
        if ismembertol(pdival,pdi_freearr,err_tol)
            pdi_index = find(pdival,pdi_freearr);
        else
            fprintf('did not find %g in pdifree within a tolerance of %g\n', pdival,err_tol)
            continue;
        end
        
        avg_across_cases(nval_index,arch_index,pdi_index) = str2double(spl_tline{avgf_col});
        nval_from_fyle(nval_index,arch_index,pdi_index) = nval;
        pdi_from_fyle(nval_index,arch_index,pdi_index) = pdival;
        num_of_cases(nval_index,arch_index,pdi_index) = num_cases(nval_index,arch_index,pdi_index) + 1;
        
        fin_vals('%g\t%d\t%s\t%s\t%g\n',pdi_from_fyle(nval_index,arch_index,pdi_index),...
            nval_from_fyle(nval_index,arch_index,pdi_index),dirstr,...
            num_of_cases(nval_index,arch_index,pdi_index)
        
        
    end % end if loop for reading lines
    
end

fclose(fin_cons);