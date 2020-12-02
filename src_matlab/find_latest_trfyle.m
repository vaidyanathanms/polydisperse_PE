function ind_latest = find_latest_trfyle(inplist,nfyles,posit) % called from output_latestfiles.m

allstrings = strsplit(char(inplist{1,1}),{'_','.'});
if isnan(str2double(allstrings{posit}))
    fprintf('Not a number at %d\t%s\n',posit,allstrings{posit});
    ind_latest = -1;
    return
end
tstamp_new = str2double(allstrings{posit});
ind_latest = 1;
for i = 2:nfyles
    allstrings = strsplit(char(inplist{i,1}),{'_','.'});
    tstamp = str2double(allstrings{posit});
    if tstamp > tstamp_new
        tstamp_new = tstamp;
        ind_latest = i;
    end
end