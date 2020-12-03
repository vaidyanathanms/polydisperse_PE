function outlist = renumber_files(inplist,numfyles) % called from adsfrac.m
outlist = cell(numfyles);
tsteparr = zeros(numfyles,1);
for i = 1:numfyles
    spl_tline = strtrim(strsplit(strtrim(inplist(i).name),{'_','.'}));
    lenline = length(spl_tline);
    tsteparr(i,1) = str2double(spl_tline(lenline-1));
end

% sort timearr
sorttstep = sort(tsteparr(:,1));

for j = 1:length(sorttstep)
    tval = sorttstep(j,1);
    for i = 1:numfyles
        spl_tline = strtrim(strsplit(strtrim(inplist(i).name),{'_','.'}));
        lenline = length(spl_tline);
        tfyle = str2double(spl_tline(lenline-1));
        if tval == tfyle
            outlist{j} = inplist(i).name;
            break;
        end
    end
end