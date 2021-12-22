%% Compute net charge within the brush

clear
clc
close all
format long

%% Color codes
green = [0 0.5 0.0]; gold = [0.9 0.75 0]; orange = [0.91 0.41 0.17];brown = [0.2 0 0];
pclr = {'r',green,'b',brown,'k','m', gold};
lsty = {'-','--',':'};
msty = {'d','o','s','x'};

%% Input data
nfreearr = [16;32;64;128;150];
casearr  = [1,2,3,4];
pdi_freearr = [1,1.5];
arch_arr = {'bl_bl','al_al'};
leg_arr  = {'Block-Block','Alter-Alter'}; % ALWAYS CHECK for correspondence with arch_arr
pdigraft = 1.0;
nmonfree = 30; nmongraft = 30; ngraft = 32;
cutoff = '1.50';
nch_graft = 32;
lz = 120; area=35^2;
set_tmax = 3e7; % maximum timestep for analysis;

s1 = create_output_dirs('./../../net_charge');

for pdi_cntr = 1:length(pdi_freearr) % begin pdi free loop
    
    ref_pdifree     = pdi_freearr(pdi_cntr);
    pdifree_str     = num2str(ref_pdifree,'%1.1f');
    
    fnr   = fopen(sprintf('./../../net_charge/QnetBound_%s.txt',pdifree_str),'w'); % case wise net charge
    fprintf(fnr,'%s \n','NetCharge: Q_{b}=\Delta(n_g)*\int(\sum(q_j n_j(z)dz, j=all entities)z=0,Lz))');
    fprintf(fnr,'%s\t%s\t%s\t%s\t%s\t%s\n','N_pa','Arch','Case','Analysis File','ht_{br}','Q_{b}');
    fmon  = fopen(sprintf('./../../net_charge/nmonlist_%s.txt',pdifree_str),'w'); % monomer distribution details
    fprintf(fmon,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n',...
        'N_pa','Arch','Case','N_gr','N_bb','N_ntbr','N_chbr','N_ntfr','N_chfr','Nposions','Nnegions','Not');
    fave  = fopen(sprintf('./../../net_charge/AvgQnetBound_%s.txt',pdifree_str),'w'); % average net charge
    fprintf(fave,'%s\t%s\t%s\t%s\t%s\n','N_pa','Arch','num_cases','<Q_{b}>','StdErrMean');
    
    pdigraft_str = num2str(pdigraft,'%1.1f');
    
    for ncnt = 1:length(nfreearr) % begin nfree loop
        nval = nfreearr(ncnt);
        
        for arch_cnt = 1:length(arch_arr)  % begin arch loop
            dirstr = arch_arr{arch_cnt};
            sumq_arr = zeros(1,1);
            sumq_across_cases = 0; q_ncases = 0;
            
            for casecntr = 1:length(casearr) % begin case loop
                casenum = casearr(casecntr);
                
                fprintf('Analyzing pdi/nfree/arch/case: %g\t%d\t%s\t%d\n',ref_pdifree,nval,dirstr,casenum);
                % Check input file name to count number of monomers of
                % each type
                inp_fylename = sprintf('./../../data_all_dir/n_%d/%s/pdifree%s_pdigraft_%s/Case_%d/PEinitdata.txt',...
                    nval,dirstr,pdifree_str,pdigraft_str,casenum);
                if exist(inp_fylename,'file') ~= 2
                    fprintf('ERROR: %s does not exist/empty file\n',inp_fylename);
                    continue;
                end
                idlist = [1;2;3;4;5;6;7;8]; % ids[grmon,bb_brush,neut_brush,charg_brush,neut_free,charg_free,posions,negions]...
                outmonlist = find_all_mons_name(inp_fylename,idlist);
                fprintf(fmon,'%d\t%s\t%d\t',nval,dirstr,casenum);
                for u = 1:length(outmonlist)
                    fprintf(fmon,'%d\t',outmonlist(u,1));
                end
                fprintf(fmon,'%d\n',sum(outmonlist(:,1)));
                
                % Now analyze all the charge details -- analyze the
                % latest time stamp
                dirname = sprintf('./../../sim_results/outresults_dir_n_%d/%s/pdifree_%s_pdigraft_%s/Case_%d',...
                    nval,dirstr,pdifree_str,pdigraft_str,casenum);
                if ~exist(dirname,'dir')
                    fprintf('%s does not exist\n',dirname);
                    continue
                end
                qnet_prefix = 'dens_config_*.lammpstrj';
                qnet_fylelist = dir(strcat(dirname,'/',qnet_prefix));
                if min(size(qnet_fylelist)) == 0
                    fprintf('No files/Empty files are found for %s\n',qnet_prefix);
                    continue;
                end
                
                nfyles = numel(qnet_fylelist); %number of files of the type
                if nfyles == 0
                    fprintf('Did not find files of the type dens_config* in %d\t%s\n', nfree,dirstr);
                    continue;
                else
                    [latest_fyleindex] = find_latest_fyle(qnet_fylelist,nfyles);
                end
                
                qnet_fylename = strcat(dirname,'/',qnet_fylelist(latest_fyleindex).name);
                if exist(qnet_fylename,'file') ~= 2
                    fprintf('%s does not exist/empty file\n',qnet_fylename);
                    continue;
                elseif struct(dir(qnet_fylename)).bytes == 0
                    fprintf('Empty file: %s \n',qnet_fylename);
                    continue;
                end
                fprintf('Analyzing netcharge using %s\n', qnet_fylename);
                
                all_data = importdata(qnet_fylename,' ',1);
                fld = all_data.data;
                rdata = fld(:,1); lz = fld(1,1) + fld(length(fld(:,1)),1); %because of binning
                
                nbins = length(rdata(:,1));
                chargearr = [0;0;1;0;-1;1;-1]; %[bb_br,nt_br,ch_br,nt_fr,ch_fr,posi,negi]
                ntotarr = sum(outmonlist); %graft has zero charge. outmonlist(1,1) is the graft
                pbackbon_brush = fld(:,2)*outmonlist(2,1)*lz*area/nbins;
                pneutral_brush = fld(:,3)*outmonlist(3,1)*lz*area/nbins;
                pcharg_brush   = fld(:,4)*outmonlist(4,1)*lz*area/nbins;
                pneutral_free  = fld(:,5)*outmonlist(5,1)*lz*area/nbins;
                pcharg_free    = fld(:,6)*outmonlist(6,1)*lz*area/nbins;
                pposions       = fld(:,7)*outmonlist(7,1)*lz*area/nbins;
                pnegions       = fld(:,8)*outmonlist(8,1)*lz*area/nbins;
                
                netcharge = chargearr(1)*pbackbon_brush + chargearr(2)*pneutral_brush ...
                    + chargearr(3)*pcharg_brush + chargearr(4)*pneutral_free + ...
                    chargearr(5)*pcharg_free + chargearr(6)*pposions + chargearr(7)*pnegions;
                
                sumq = 0.5*netcharge(1);
                
                % Find height of the brush to compute charge within the
                % brush region
                grp_prefix = 'polydens_config_*.lammpstrj';
                grp_fylelist = dir(strcat(dirname,'/',grp_prefix));
                if min(size(grp_fylelist)) == 0
                    fprintf('No files/Empty files are found for %s\n',grp_prefix);
                    continue;
                end
                nfyles = numel(grp_fylelist); %number of files of the type
                if nfyles == 0
                    fprintf('Did not find files of the type dens_config* in %d\t%s\n', nval,dirstr);
                    continue
                else
                    [latest_fyleindex] = find_latest_fyle(grp_fylelist,nfyles);
                end
                grp_fylename = strcat(dirname,'/',grp_fylelist(latest_fyleindex).name);
                [i_edge, brht] = find_brush_height(grp_fylename,rdata);
                if i_edge == -1
                    fprintf('ERROR: Did not find edge of the brush in %d\t%s\t%d\n',nval,dirstr,casenum);
                    continue
                end
                qofr = zeros(i_edge,1); rqofr = zeros(i_edge,1);
                qofr(1,1) = sumq; rqofr(1,1) = 0.5*rdata(1);
                
                % Integrate charge to the edge of the brush
                for k = 2:i_edge
                    
                    sumq = sumq + 0.5*(rdata(k)-rdata(k-1))*(netcharge(k)+netcharge(k-1));
                    qofr(k,1)  = sumq;
                    rqofr(k,1) = 0.5*(rdata(k)+rdata(k+1));
                    
                end
                
                fprintf(fnr,'%d\t%s\t%d\t%s\t%g\t%g\n',nval,dirstr,casenum,grp_fylelist(latest_fyleindex).name,brht,sumq);
                sumq_across_cases = sumq_across_cases + sumq;
                q_ncases = q_ncases + 1;
                sumq_arr(casecntr,1) = sumq;
                
            end % end case loop
            
            sumq_nozero = sumq_arr(sumq_arr(:,1)~=0);
            if q_ncases ~= length(sumq_nozero)
                fprintf('ERR: Not equal number of cases from two methods: %d\t%d\n', q_ncases,length(sumq_nozero));
                fprintf(fave,'%d\t%s\t%s\t%s\t%s\n',nval,dirstr,'ERR','ERR','ERR');
                continue
            end
            stderr = std(sumq_nozero)/sqrt(length(sumq_nozero));
            
            fprintf(fave,'%d\t%s\t%d\t%g\t%g\n',nval,dirstr,q_ncases,sumq_across_cases/q_ncases,stderr);
            
        end % end arch loop
        
    end % end ncnt loop
    
    fclose(fave); fclose(fnr); fclose(fmon);
    
end % end pdi loop

