function [] = make_grid_winds_zonal(go_latm,go_late,go_mask,str_nameout,par_tauopt,par_plotformat)
% make_grid_winds_zonal
%
%   *********************************************************
%   *** CREATE ZONAL AVERAGE WINDSTRESS                   ***
%   *********************************************************
%
%   make_windstress(DUM_NI,DUM_NJ,loc_nswl)
%   creates a synthetic zonal average windstress and takes 3 arguments:
%
%   DUM_NI [INTEGER] (e.g. 36)
%   --> the i-dimension of the grid
%   DUM_NJ [INTEGER] (e.g. 36)
%   --> the j-dimension of the grid
%   loc_nswl [STRING] (e.g. 'nwsw')
%   --> the profile required in N and S hemispheres
%       valid options are:
%       'nwsw' == NH waterworld, SH waterworld
%       'nlsw' == NH land present, SH waterworld
%       'nwsl' == NH waterworld, SH land present
%       'nlsl' == NH land present, SH land present
%
%   *********************************************************
%
%   ***********************************************************************
%   *** HISTORY ***********************************************************
%   ***********************************************************************
%
%   14/09/25: CREATED
%   14/09/26: initial working version completed
%   14/10/07: extended for get tau on both u and v grid
%   14/10/28: added wind
%   17/02/23: adapted for inforporation muffingen & renamed
%   17/04/18: revised automatic detection of gateways
%             revised new algorithm for calculating wind stress
%             removed un-used input parameters
%   17/10/22: added new input parameter
%   18/09/11: added minimum wind stress
%   20/02/26: added 'grey world' intermediate stress option (3)
%             re-ordered land vs. ocean array index to be consistent 
%             with the actual choice options (e.g. 1 == land)
%   24/03/13: added windspeed output for usage with ENTS
%
%   ***********************************************************************

% *********************************************************************** %
% *** INITIALIZE PARAMETERS & VARIABLES ********************************* %
% *********************************************************************** %
%
% START
% determine grid size
[jmax imax] = size(go_mask);
switch par_tauopt
    case {1}
        % land world (modern NH / paleo Eocene (both hemispheres))
        loc_nswl(1:2) = 'nl';
        loc_nswl(3:4) = 'sl';
        % report action
        disp(['         -> Zonal wind stress profile to be generated consistent with land barriers in both hemisphere.']);
    case {2}
        % water world
        loc_nswl(1:2) = 'nw';
        loc_nswl(3:4) = 'sw';
        % report action
        disp(['         -> Zonal wind stress profile to be generated consistent with no barriers in eithher hemisphere.']);
    case {3}
        % grey world
        loc_nswl(1:2) = 'ng';
        loc_nswl(3:4) = 'sg';
        % report action
        disp(['         -> An intermediate (between having land and no land barriers) strength zonal wind stress profile will be generated.']);
    otherwise
        % automatically determine water/land options
        % NOTE: any zeros present in mask == non-ocean => land present
        %       <=> a '1.0' will indicate a zonal gateway in the mask mean
        % NOTE: assume a higher lat than 40 degree N/S criteria
        % NOTE: remember that '1' is ocean
        % NOTE: re-orientate the latitude vector (so N pole is at the 'top')
        vo_latm = fliplr(go_latm);
        lat_thrsh = 40.0;
        vo_mask  = mean(go_mask');
        % report action
        disp(['         -> cookiegen will very crudely try and identify the presence of high latitude ocean gateways and set a corresponding wind stress profile.']);
        % N
        if max(vo_mask(find(vo_latm > lat_thrsh))) < 1.0
            loc_nswl(1:2) = 'nl';
            fprintf('       * No Northern gateway found.\n')
        else
            loc_nswl(1:2) = 'nw';
            fprintf('       * Northern gateway (> 40N) found.\n')
        end
        % S
        if max(vo_mask(find(vo_latm < -lat_thrsh))) < 1.0
            loc_nswl(3:4) = 'sl';
            fprintf('       * No Southern gateway found.\n')
        else
            loc_nswl(3:4) = 'sw';
            fprintf('       * Southern gateway (< -40S) found.\n')
        end
end
% set GOLDSTEIn parameters :: drag coefficent
go_cd = 0.0013;
% set GOLDSTEIn parameters :: air density
go_rhoair = 1.25;
% set zonal profile parameters -- with land
% NOTE: derived from NH observations
%       (and to be vaguely consistent with HadCM3 Eocene simulations)
par_ws_amp(1)         = -0.090;
par_ws_offset(1)      = -0.020;
par_ws_cyclefactor(1) = 6.0;
par_ws_modfactor(1)   = 2.0;
par_ws_modpower(1)    = 2.0;
% set zonal profile parameters -- water world
% NOTE: derived vaguely from 
%       Marshall et al. [2007] and Enderton and Marshall [2014]
%       also some consistency with SH observations
par_ws_amp(2)         = -0.140;
par_ws_offset(2)      = -0.020;
par_ws_cyclefactor(2) = 6.0;
par_ws_modfactor(2)   = 2.0;
par_ws_modpower(2)    = 2.0;
% set zonal profile parameters -- grey (intermediate)
par_ws_amp(3)         = -0.115;
par_ws_offset(3)      = -0.020;
par_ws_cyclefactor(3) = 6.0;
par_ws_modfactor(3)   = 2.0;
par_ws_modpower(3)    = 2.0;
% create zonal output arrays
taux_u = zeros(jmax,imax);
taux_u_1d = zeros(jmax,1);
taux_v = zeros(jmax,imax);
taux_v_1d = zeros(jmax,1);
% set minimum tau (N m-2)
% NOTE: derived from JOSEY et al. [2002]
par_tau_min = 0.02;
% set date
str_date = [datestr(date,11), datestr(date,5), datestr(date,7)];
%
% *********************************************************************** %

% *********************************************************************** %
% *** CREATE ZONAL WIND-STRESS PROFILE ********************************** %
% *********************************************************************** %
%
% NOTE: TAKE CARE WITH C-GRID!!!
%
%   ----v----
%   |       |
%   |   r   u
%   |       |
%   ---------
%
% Southern hemisphere
for j=1:jmax/2
    switch loc_nswl
        case {'nwsl', 'nlsl'}
            loc_amp         = par_ws_amp(1);
            loc_offset      = par_ws_offset(1);
            loc_cyclefactor = par_ws_cyclefactor(1);
            loc_modfactor   = par_ws_modfactor(1);
            loc_modpower    = par_ws_modpower(1);
        case {'nwsw', 'nlsw'}
            loc_amp         = par_ws_amp(2);
            loc_offset      = par_ws_offset(2);
            loc_cyclefactor = par_ws_cyclefactor(2);
            loc_modfactor   = par_ws_modfactor(2);
            loc_modpower    = par_ws_modpower(2);
        otherwise
            loc_amp         = par_ws_amp(3);
            loc_offset      = par_ws_offset(3);
            loc_cyclefactor = par_ws_cyclefactor(3);
            loc_modfactor   = par_ws_modfactor(3);
            loc_modpower    = par_ws_modpower(3);
    end
    % SIN((3.1416/180)*$L$5*ABS(C10))
    loc_tmp1  = sin( loc_cyclefactor*(pi/180.0)*abs(go_latm(j)) );
    loc_tmp1e = sin( loc_cyclefactor*(pi/180.0)*abs(go_late(j+1)) );
    % SIN((3.1416/180)*$L$3*ABS(C10))
    loc_tmp2  = sin( loc_modfactor*(pi/180.0)*abs(go_latm(j)) );
    loc_tmp2e = sin( loc_modfactor*(pi/180.0)*abs(go_late(j+1)) );
    % $L$2*(ABS($F22))^$L$6*$E22 + $L$4*(90-ABS($C22))/90
    taux_u(j,:) = loc_amp*loc_tmp1*loc_tmp2^loc_modpower   + loc_offset*(90.0-abs(go_latm(j)))/90.0;
    taux_v(j,:) = loc_amp*loc_tmp1e*loc_tmp2e^loc_modpower + loc_offset*(90.0-abs(go_late(j+1)))/90.0;
end
% Northern hemisphere
for j=(jmax/2+1):jmax
    switch loc_nswl
        case {'nlsw', 'nlsl'}
            loc_amp         = par_ws_amp(1);
            loc_offset      = par_ws_offset(1);
            loc_cyclefactor = par_ws_cyclefactor(1);
            loc_modfactor   = par_ws_modfactor(1);
            loc_modpower    = par_ws_modpower(1);
        case {'nwsw', 'nwsl'}
            loc_amp         = par_ws_amp(2);
            loc_offset      = par_ws_offset(2);
            loc_cyclefactor = par_ws_cyclefactor(2);
            loc_modfactor   = par_ws_modfactor(2);
            loc_modpower    = par_ws_modpower(2);
        otherwise
            loc_amp         = par_ws_amp(3);
            loc_offset      = par_ws_offset(3);
            loc_cyclefactor = par_ws_cyclefactor(3);
            loc_modfactor   = par_ws_modfactor(3);
            loc_modpower    = par_ws_modpower(3);
    end
    % SIN((3.1416/180)*$L$5*ABS(C10))
    loc_tmp1  = sin( loc_cyclefactor*(pi/180.0)*abs(go_latm(j)) );
    loc_tmp1e = sin( loc_cyclefactor*(pi/180.0)*abs(go_late(j+1)) );
    % SIN((3.1416/180)*$L$3*ABS(C10))
    loc_tmp2  = sin( loc_modfactor*(pi/180.0)*abs(go_latm(j)) );
    loc_tmp2e = sin( loc_modfactor*(pi/180.0)*abs(go_late(j+1)) );
    % $L$2*(ABS($F22))^$L$6*$E22 + $L$4*(90-ABS($C22))/90
    taux_u(j,:) = loc_amp*loc_tmp1*loc_tmp2^loc_modpower   + loc_offset*(90.0-abs(go_latm(j)))/90.0;
    taux_v(j,:) = loc_amp*loc_tmp1e*loc_tmp2e^loc_modpower + loc_offset*(90.0-abs(go_late(j+1)))/90.0;
end
% impose minimum magnitude of wind stress
% NOTE: default value of par_tau_min from JOSEY et al. [2002]
loc_ans = intersect(find(abs(taux_u(:,:)) < par_tau_min),find(taux_u(:,:) > 0.0));
taux_u(loc_ans) = par_tau_min;
loc_ans = intersect(find(abs(taux_u(:,:)) < par_tau_min),find(taux_u(:,:) < 0.0));
taux_u(loc_ans) = -par_tau_min;
loc_ans = intersect(find(abs(taux_v(:,:)) < par_tau_min),find(taux_v(:,:) > 0.0));
taux_v(loc_ans) = par_tau_min;
loc_ans = intersect(find(abs(taux_v(:,:)) < par_tau_min),find(taux_v(:,:) < 0.0));
taux_v(loc_ans) = -par_tau_min;
%
% *********************************************************************** %

% *********************************************************************** %
% *** CREATE WINDS ****************************************************** %
% *********************************************************************** %
%
% derive wind speeds on t-grid (i.e. the v grid point for x)
% NOTE: from BIOGEM code :: fun_calc_u(i,j) = sqrt((sqrt(tv**2 + tv2**2))/(goldstein_rhoair*goldstein_cd))
wind_u = sign(taux_u(:,:)).*((taux_u(:,:).^2).^0.5/(go_rhoair*go_cd)).^0.5;
%
% *** CREATE 2D WINDS OCEAN vs LAND ************************************* %
%
% Windspeed over land tends to be (~2x) weaker than windspeed over ocean
% NOTE: With ENTS, need for 2D ocean-land defined windspeed
%
% make land grid (land=1, ocn=0)
go_maskl = abs(go_mask-1);
% make winds over land 2x as weak as over ocean
wind_u_ents = go_maskl.*(wind_u./2);
% fill remainder with ocean wind speed
wind_u_ents(wind_u_ents==0) = wind_u(wind_u_ents==0);
% windspeed is non-directional. remove sign
wind_u_ents = abs(wind_u_ents);
%
% *********************************************************************** %

% *********************************************************************** %
% *** OUTPUT RESULTS **************************************************** %
% *********************************************************************** %
%
% *** PLOT PROFILES ***************************************************** %
%
% calculate zonal mean (should be equal to any particular longitude!)
for j=1:jmax
    taux_u_1d(j) = mean(taux_u(j,:));
    taux_v_1d(j) = mean(taux_v(j,:));
end
% plot figure -- u grid
if (isempty(par_plotformat))
    %
else
    figure;
    plot(go_latm(:),taux_u_1d(:));
    axis([-90 90 -0.1 0.2]);
    hold on;
    scatter(go_latm(:),taux_u_1d(:));
    % save figure
    if (strcmp(par_plotformat,'pdf'))
        exportgraphics(gcf,[str_nameout, '.taux_u.' str_date '.pdf'],'BackgroundColor','none','ContentType','vector');
    else
        exportgraphics(gcf,[str_nameout, '.taux_u.' str_date '.' par_plotformat]);
    end
end
% plot figure -- v grid
if (isempty(par_plotformat))
    %
else
    figure;
    plot(go_late(2:end),taux_v_1d(:));
    axis([-90 90 -0.1 0.2]);
    hold on;
    scatter(go_late(2:end),taux_v_1d(:));
    % save figure
    if (strcmp(par_plotformat,'pdf'))
        exportgraphics(gcf,[str_nameout, '.taux_v.' str_date '.pdf'],'BackgroundColor','none','ContentType','vector');
    else
        exportgraphics(gcf,[str_nameout, '.taux_v.' str_date '.' par_plotformat]);
    end
end
%
% *** SAVE FILES ******************************************************** %
%
% open file
fid = fopen([str_nameout, '.taux_u.dat'],'w');
% write data
for j=1:jmax
    for i=1:imax
        fprintf(fid,'  %d',taux_u(j,i));
        fprintf(fid,'\n');
    end
end
% close file
fclose(fid);
% open file
fid = fopen([str_nameout, '.tauy_u.dat'],'w');
% write data
for j=1:jmax
    for i=1:imax
        fprintf(fid,'  %d',0.0);
        fprintf(fid,'\n');
    end
end
% close file
fclose(fid);
%
% open file
fid = fopen([str_nameout, '.taux_v.dat'],'w');
% write data
for j=1:jmax
    for i=1:imax
        fprintf(fid,'  %d',taux_v(j,i));
        fprintf(fid,'\n');
    end
end
% close file
fclose(fid);
% open file
fid = fopen([str_nameout, '.tauy_v.dat'],'w');
% write data
for j=1:jmax
    for i=1:imax
        fprintf(fid,'  %d',0.0);
        fprintf(fid,'\n');
    end
end
% close file
fclose(fid);
%
% open file
fid = fopen([str_nameout, '.wvelx.dat'],'w');
% write data
for j=1:jmax
    for i=1:imax
        fprintf(fid,'  %d',wind_u(j,i));
        fprintf(fid,'\n');
    end
end
% close file
fclose(fid);
% open file
fid = fopen([str_nameout, '.wvely.dat'],'w');
% write data
for j=1:jmax
    for i=1:imax
        fprintf(fid,'  %d',0.0);
        fprintf(fid,'\n');
    end
end
% close file
fclose(fid);
%
% Save 2D ASCII windspeed scalar on full grid
outname = [str_nameout, '.windspeed__ents.dat'];
a = wind_u_ents;
save(outname,'a','-ascii');
%
% *********************************************************************** %

% *********************************************************************** %
% *** END *************************************************************** %
% *********************************************************************** %
%
% END
%
% *********************************************************************** %
