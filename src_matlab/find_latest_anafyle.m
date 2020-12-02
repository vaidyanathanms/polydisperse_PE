function ind_latest = find_latest_anafyle(inplist,nfyles,posit) % called from output_latestfiles.m

allstrings = strsplit(inplist(1).name,{'_','.'});
if isnan(str2double(allstrings{posit}))
    fprintf('Not a number at %d\t%s\n',posit,allstrings{posit});
    ind_latest = -1;
    return
end
tstamp_new = str2double(allstrings{posit});
ind_latest = 1;
for i = 2:nfyles
    allstrings = strsplit(inplist(i).name,{'_','.'});
    tstamp = str2double(allstrings{posit});
    if tstamp > tstamp_new
        tstamp_new = tstamp;
        ind_latest = i;
    end
end