function ind_latest = find_latest_fyle(inplist,nfyles)

fylename = strsplit(inplist(1).name,{'_','.'});
tstamp_new = str2double(fylename{3});
ind_latest = 1;
for i = 1:nfyles
    fylename = strsplit(inplist(i).name,{'_','.'});
    tstamp = str2double(fylename{3});
    if tstamp > tstamp_new
        tstamp_new = tstamp;
        ind_latest = i;
    end
end

