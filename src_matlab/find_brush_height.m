function [bhtindex,bht] = find_brush_height(grpfyle,dens_rdata)
bhtindex = -1; bht = 0;
if exist(grpfyle,'file') ~= 2
    fprintf('%s does not exist/empty file\n',grpfyle);
    return;
elseif struct(dir(grpfyle)).bytes == 0
    fprintf('Empty file: %s \n',qnet_fylename);
    return;
end

fid_g  = fopen(grpfyle);
data_g = textscan(fid_g,'%f%f%f','Headerlines',1);

fld_g   = cell2mat(data_g);
rgrp    = fld_g(:,1);
dens_g  = fld_g(:,2);
[maxden, imaxden]  = max(dens_g);

% Find edge of brush
for k = imaxden:length(dens_g)
    if dens_g(k,1) < 0.05*maxden
        r_edge = rgrp(k,1);
        break;
    end
end

% Find index corresponding to r_edge in dens_rdata
for u = 1:length(dens_rdata)
    if dens_rdata(u,1) == r_edge
        bhtindex = u; bht = dens_rdata(u,1);
        break;
    elseif dens_rdata(u,1) > r_edge
        bhtindex = u; bht = 0.5*(dens_rdata(u,1)+dens_rdata(u-1,1));
        break;
    end
end
