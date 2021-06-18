function timeout = extract_timeval(inpstr)
splstr = strsplit(inpstr,{'_','.'});
lenstr = length(splstr);
timeout = str2double(splstr{lenstr-1});