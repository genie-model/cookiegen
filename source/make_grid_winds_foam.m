function [] = make_grid_winds_foam(gilonpm,gilonpe,gilatpm,gilatpe,gilonam,gilonae,gilatam,gilatae,gimask,golonm,golone,golatm,golate,gomask,str,par_plotformat)
% make_grid_winds_foam
%
%   *********************************************************
%   *** RE_GRID WIND PRODUCTS                             ***
%   *********************************************************
%
%
%   wstress is organised :   level 1   tau_x at u point
%                            level 2   tau_x at v point
%                            level 3   tau_y at u point
%                            level 4   tau_y at v point
%
%   wvelocity is organised : level 1   x velocity at grid point
%                            level 2   y velocity at grid point
%
%   wspeed is organised :    level 1   speed at grid point
%
%   str input KEY:
%   str(1).nc == par_nc_axes_name
%   str(2).nc == par_nc_topo_name
%   str(3).nc == par_nc_atmos_name
%   str(4).nc == par_nc_ocean_name
%   str(5).nc == par_nc_coupl_name
%
%   *********************************************************
%
%   ***********************************************************************
%%

% *********************************************************************** %
% *** INITIALIZE ******************************************************** %
% *********************************************************************** %
%
% determine output grid size (remember: [rows columns])
[jmax, imax] = size(gomask);
% simple parsing of par_plotformat compatible with existing cookiegen code
if (isempty(par_plotformat))
    optplots = false;
elseif (strcmp(par_plotformat,'pdf'))
    optplots = true;
else
    optplots = false;
end
%
% *********************************************************************** %

% *********************************************************************** %
% *** SET UP GRID ******************************************************* %
% *********************************************************************** %
%
% NOTE: the u, v winds are on a different grid to the wind stress
%
% *** load data -- wind stress ****************************************** %
%
% open netCDF file
ncid = netcdf.open([str(1).path '/' str(1).exp '/' str(5).nc '.nc'],'nowrite');
% read netCDf information
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
% load TAUX
varid  = netcdf.inqVarID(ncid,'TAUX');
loctaux(:,:,:) = netcdf.getVar(ncid,varid);
% load TAUY
varid  = netcdf.inqVarID(ncid,'TAUY');
loctauy(:,:,:) = netcdf.getVar(ncid,varid);
% close netCDF file
netcdf.close(ncid);
%
% *** process data -- wind stress *************************************** %
%
% create annual averages
giwtauu = 0.0*loctaux(:,:,1);
giwtauv = 0.0*loctauy(:,:,1);
for t=1:12,
    giwtauu = giwtauu + loctaux(:,:,t)/12.0;
    giwtauv = giwtauv + loctauy(:,:,t)/12.0;
end
% re-orientate
giwtauu = flipud(giwtauu');
giwtauv = flipud(giwtauv');
%
% *** load data -- wind velocity **************************************** %
%
% open netCDF file
ncid = netcdf.open([str(1).path '/' str(1).exp '/' str(3).nc '.nc'],'nowrite');
% read netCDf information
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
% load U
varid  = netcdf.inqVarID(ncid,'U');
locwvelu(:,:,:,:) = netcdf.getVar(ncid,varid);
% load V
varid  = netcdf.inqVarID(ncid,'V');
locwvelv(:,:,:,:) = netcdf.getVar(ncid,varid);
% close netCDF file
netcdf.close(ncid);
%
% *** process data -- wind velocity ************************************* %
%
% NOTE: atmosphere has 18 levels -- #18 is 'surface' (992 mbar)
% create annual averages
giwvelu = 0.0*locwvelu(:,:,18,1);
giwvelv = 0.0*locwvelv(:,:,18,1);
for t=1:12
    giwvelu = giwvelu + locwvelu(:,:,18,t)/12.0;
    giwvelv = giwvelv + locwvelv(:,:,18,t)/12.0;
end
% re-orientate
giwvelu = flipud(giwvelu');
giwvelv = flipud(giwvelv');
%
% *** process data -- wind speed **************************************** %
%
% calculate wind speed -- annual average from uv
% (for backwards-compatability) 
giwspd_uvaa = (giwvelu.^2 + giwvelv.^2).^0.5;
% calculate wind speed -- monthly average from uv
% [better averaging] 
giwspd_uvma = zeros(size(locwvelu(:,:,18,1)));
for t=1:12
    giwspd_uvma = giwspd_uvma + (locwvelu(:,:,18,t).^2 + locwvelv(:,:,18,t).^2).^0.5/12.0;
end
% re-orientate
giwspd_uvma = flipud(giwspd_uvma');
%
% NOTE: windspeed directly from output (wsma) is unavailable for FOAM
%
% *** load mask -- wind velocity **************************************** %
%
% open netCDF file
ncid = netcdf.open([str(1).path '/' str(1).exp '/' str(3).nc '.nc'],'nowrite');
% read netCDf information
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
% load MASK
varid  = netcdf.inqVarID(ncid,'ORO');
locmaskvel(:,:,:) = netcdf.getVar(ncid,varid);
% close netCDF file
netcdf.close(ncid);
%
% *** process data -- wind velocity ************************************ %
%
% create annual averages (of mask!!!)
gimaskpl = 0.0*locmaskvel(:,:,1);
for t=1:12
    gimaskpl = gimaskpl + locmaskvel(:,:,t)/12.0;
end
% re-orientate
gimaskpl = flipud(gimaskpl');
% derive mask values (ocean == 1, land == NaN, seaice > 1)
% NOTE: FOAM mask is odd ... land is not quite equal to 1,
%       more like 0.999999999 something ... unclear!
% create land (only) mask (land==1, ocean/seaice==NaN)
gimaskpl(find(gimaskpl<0.999)) = NaN;
gimaskpl(find(gimaskpl>1.0))   = NaN;
gimaskpl(find(~isnan(gimaskpl))) = 1.0;
%
% *** load mask -- wind stress ****************************************** %
%
% open netCDF file
ncid = netcdf.open([str(1).path '/' str(1).exp '/' str(5).nc '.nc'],'nowrite');
% read netCDf information
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
% load MASK
varid  = netcdf.inqVarID(ncid,'ORO');
locmask(:,:,:) = netcdf.getVar(ncid,varid);
% close netCDF file
netcdf.close(ncid);
%
% *** process data -- wind stress mask ********************************** %
%
% create annual averages (of mask!!!)
gimaskp = 0.0*locmask(:,:,1);
for t=1:12,
    gimaskp = gimaskp + locmask(:,:,t)/12.0;
end
% re-orientate
gimaskp = flipud(gimaskp');
% derive mask values (ocean == 1, land == NaN, seaice > 1)
% NOTE: FOAM mask is odd ... land is not quite equal to 1,
%       more like 0.999999999 something ... unclear!
gimaskp(find(gimaskp<0.999)) = 1.0;
gimaskp(find(gimaskp>1.0))   = 1.0;
gimaskp(find(gimaskp~=1.0))  = NaN;
%
% *********************************************************************** %

% *********************************************************************** %
% *** RE-GRID *********************************************************** %
% *********************************************************************** %
%
% *** set up GOLDSTEIN grid ********************************************* %
%
% Need to determine wind stress at u and v points of the Arakawa C
% grid in GENIE ...
%
%   ----v----
%   |       |
%   |   c   u
%   |       |
%   ---------
% 
% remember: [rows columns] == [j i]
% grid doudaries are as follows:
% GENIE c-grid: (golate,golone)   == jmax+1 x imax+1
% GENIE u-grid: (golatue,golonue) == jmax   x imax+1
% GENIE v-grid: (golatve,golonve) == jmax+1 x imax
%
% create GENIE u and v edge axes
% latitude
golatue = golate;
golatve = [golatm 90.0];
% longitude
golonue = [golonm golonm(end)+360.0/imax];
golonve = golone;
% create GENIE NaN mask -- ocean
gm = gomask;
gm(find(gm == 0)) = NaN;
% create GENIE NaN mask -- land
gml = gomask;
gml(find(gml == 1)) = NaN;
gml(find(gml == 0)) = 1;
%
% *** Put on a GOLDSTEIN grid ******************************************* %
%
% NOTE: don't forget to flip the re-gridded orientation back around again
%
% plot raw wind stress
if (optplots), plot_2dgridded(flipud(giwtauu),999.0,'',[[str(2).dir '/' str(2).exp] '.wtau_u.IN'],['wind stress in -- u']); end
if (optplots), plot_2dgridded(flipud(giwtauv),999.0,'',[[str(2).dir '/' str(2).exp] '.wtau_v.IN'],['wind stress in -- v']); end
% 
% -> wind stress @ u point
% NOTE: GENIE u-grid: (golatue,golonue) == jmax   x imax+1
% apply p-grid mask to input wind stress grid
% replace NaNs with zeros
% apply GENIE mask to output wind stress
% u
fprintf('       - Regridding wind stress (x) to GOLDSTEIN u-grid\n');
[gowtauuu,gofwtauuu] = make_regrid_2d(gilonpe,gilatpe,(gimaskp.*giwtauu)',golonue,golatue,false);
gowtauuu(find(isnan(gowtauuu))) = 0.0;
gowtauuu = gowtauuu'; 
gowtauuu = gomask.*gowtauuu;
if (optplots), plot_2dgridded(flipud(gm.*gowtauuu),999.0,'',[[str(2).dir '/' str(2).exp] '.wtau_xATu.out'],['wind stress out -- x @ u']); end
% v
fprintf('       - Regridding wind stress (y) to GOLDSTEIN u-grid\n');
[gowtauvu,gofwtauvu] = make_regrid_2d(gilonpe,gilatpe,(gimaskp.*giwtauv)',golonue,golatue,false);
gowtauvu(find(isnan(gowtauvu))) = 0.0;
gowtauvu = gowtauvu';
gowtauvu = gomask.*gowtauvu;
if (optplots), plot_2dgridded(flipud(gm.*gowtauvu),999.0,'',[[str(2).dir '/' str(2).exp] '.wtau_yATu.out'],['wind stress out -- y @ u']); end
% 
% -> wind stress @ v point
% NOTE: GENIE v-grid: (golatve,golonve) == jmax+1 x imax
% apply p-grid mask to input wind stress grid
% replace NaNs with zeros
% apply GENIE mask to output wind stress
% u
fprintf('       - Regridding wind stress (x) to GOLDSTEIN v-grid\n');
[gowtauuv,gofwtauuv] = make_regrid_2d(gilonpe,gilatpe,(gimaskp.*giwtauu)',golonve,golatve,false);
gowtauuv(find(isnan(gowtauuv))) = 0.0;
gowtauuv = gowtauuv'; 
gowtauuv = gomask.*gowtauuv;
if (optplots), plot_2dgridded(flipud(gm.*gowtauuv),999.0,'',[[str(2).dir '/' str(2).exp] '.wtau_xATv.out'],['wind stress out -- x @ v']); end
% v
fprintf('       - Regridding wind stress (y) to GOLDSTEIN v-grid\n');
[gowtauvv,gofwtauvv] = make_regrid_2d(gilonpe,gilatpe,(gimaskp.*giwtauv)',golonve,golatve,false);
gowtauvv(find(isnan(gowtauvv))) = 0.0;
gowtauvv = gowtauvv';
gowtauvv = gomask.*gowtauvv;
if (optplots), plot_2dgridded(flipud(gm.*gowtauvv),999.0,'',[[str(2).dir '/' str(2).exp] '.wtau_yATv.out'],['wind stress out -- y @ v']); end
% 
% -> wind velocity
% NOTE: GENIE c-grid: (golate,golone) == jmax+1 x imax+1
% NOTE: DO NOT apply a mask (as the output is used globally by the EMBM)
% NOTE: input grid is the FOAM atmosphere grid
% replace NaNs with zeros
% u
fprintf('       - Regridding wind velocity (x) to GOLDSTEIN c-grid\n');
if (optplots), plot_2dgridded(flipud(giwvelu),999.0,'',[[str(2).dir '/' str(2).exp] '.wvel_x.IN'],['wind velocity in -- x']); end
[gowvelu,gofwvelu] = make_regrid_2d(gilonae,gilatae,giwvelu',golone,golate,false);
gowvelu(find(isnan(gowvelu))) = 0.0;
gowvelu = gowvelu'; 
gofwvelu = gofwvelu';
if (optplots), plot_2dgridded(flipud(gowvelu),999.0,'',[[str(2).dir '/' str(2).exp] '.wvel_x.OUT'],['wind velocity out -- x']); end
% v
fprintf('       - Regridding wind velocity (y) to GOLDSTEIN c-grid\n');
if (optplots), plot_2dgridded(flipud(giwvelv),999.0,'',[[str(2).dir '/' str(2).exp] '.wvel_y.IN'],['wind velocity in -- y']); end
[gowvelv,gofwvelv] = make_regrid_2d(gilonae,gilatae,giwvelv',golone,golate,false);
gowvelv(find(isnan(gowvelv))) = 0.0;
gowvelv = gowvelv';
gofwvelv = gofwvelv';
if (optplots), plot_2dgridded(flipud(gowvelv),999.0,'',[[str(2).dir '/' str(2).exp] '.wvel_y.OUT'],['wind velocity out -- y']); end
%
% -> wind speed
% NOTE: GENIE c-grid: (golate,golone) == jmax+1 x imax+1
% NOTE: input grid is the FOAM atmosphere grid
% NOTE: create only no-mask version
% replace NaNs with zeros
% apply GENIE mask to output wind speed
fprintf('       - Regridding wind speed to GOLDSTEIN c-grid\n');
% wspeed_uvaa
plot_2dgridded(flipud(giwspd_uvaa),999.0,'',[[str(2).dir '/' str(2).exp] '.wspd_uvaa.IN'],['wind speed in']);
[gowspdall,gofwspdall] = make_regrid_2d(gilonae,gilatae,giwspd_uvaa',golone,golate,false);
gowspdall = gowspdall';
gofwspdall = gofwspdall';
if (optplots), plot_2dgridded(flipud(gowspdall),999.0,'',[[str(2).dir '/' str(2).exp] '.wspd_uvaa.OUTALL'],['wind speed out -- no mask']); end
wspeed_uvaa = gomask.*gowspdall;
wspeed_uvaa(find(isnan(wspeed_uvaa))) = 0.0;
% wspeed_uvma
plot_2dgridded(flipud(giwspd_uvma),999.0,'',[[str(2).dir '/' str(2).exp] '.wspd_uvma.IN'],['wind speed in']);
[gowspdall_uvma,gofwspdall_uvma] = make_regrid_2d(gilonae,gilatae,giwspd_uvma',golone,golate,false);
gowspdall_uvma = gowspdall_uvma';
gofwspdall_uvma = gofwspdall_uvma';
if (optplots), plot_2dgridded(flipud(gowspdall_uvma),999.0,'',[[str(2).dir '/' str(2).exp] '.wspd_uvma.OUTALL'],['wind speed out -- no mask']); end
wspeed_uvma = gomask.*gowspdall_uvma;
wspeed_uvma(find(isnan(wspeed_uvma))) = 0.0;
%
%
% *** Put on a ENTS (land) grid ******************************************* %
%
% NOTE: don't forget to flip the re-gridded orientation back around again
%
% make GENIE land mask
glmask = abs(gomask-1);
%
% -> wind speed
% NOTE: GENIE c-grid: (golate,golone) == jmax+1 x imax+1
% NOTE: input grid is the FOAM atmosphere grid
% NOTE: create only no-mask version
% replace NaNs with zeros
% apply GENIE mask to output wind speed
fprintf('       - Regridding wind speed to land c-grid\n');
% re-grid masked field wspeed_uvaa
[gowspdl,gofwspdl] = make_regrid_2d(gilonae,gilatae,(gimaskpl.*giwspd_uvaa)',golone,golate,false);
gowspdl = gowspdl';
gofwspdl = gofwspdl';
if (optplots), plot_2dgridded(flipud(gowspdl),999.0,'',[[str(2).dir '/' str(2).exp] '.wspd_uvaal.OUT'],['wind speed out -- land mask']); end
wspeed_uvaal = glmask.*gowspdl;
wspeed_uvaal(find(isnan(wspeed_uvaal))) = 0.0;
% re-grid masked field wspeed_uvma
[gowspdl,gofwspdl] = make_regrid_2d(gilonae,gilatae,(gimaskpl.*giwspd_uvma)',golone,golate,false);
gowspdl = gowspdl';
gofwspdl = gofwspdl';
if (optplots), plot_2dgridded(flipud(gowspdl),999.0,'',[[str(2).dir '/' str(2).exp] '.wspd_uvmal.OUT'],['wind speed out -- land mask']); end
wspeed_uvmal = glmask.*gowspdl;
wspeed_uvmal(find(isnan(wspeed_uvmal))) = 0.0;
%
% *** Copy to output arrays ********************************************* %
%
% -> wind stress
wstress(:,:,1) = flipud(gowtauuu); %g_taux_u
wstress(:,:,2) = flipud(gowtauuv); %g_taux_v
wstress(:,:,3) = flipud(gowtauvu); %g_tauy_u
wstress(:,:,4) = flipud(gowtauvv); %g_tauy_v
% -> wind velocity
wvelocity(:,:,1) = flipud(gowvelu);
wvelocity(:,:,2) = flipud(gowvelv);
%
% *********************************************************************** %

% *********************************************************************** %
% *** SAVE DATA ********************************************************* %
% *********************************************************************** %
%
% Save regridded data to file OCEAN-only
% Taux at u point (g_taux_u == gowtauuu)
outname = [str(2).dir '/' str(2).exp '.taux_u.dat'];
c = wstress(:,:,1); b = permute(c, [2 1]); a = reshape(b, [imax*jmax 1]);
save (outname, 'a', '-ascii');
fprintf('       - Written tau u (u point) data to %s\n',outname);
% Taux at v point (g_taux_v == gowtauuv)
outname = [str(2).dir '/' str(2).exp '.taux_v.dat'];
c = wstress(:,:,2); b = permute(c, [2 1]); a = reshape(b, [imax*jmax 1]);
save (outname, 'a', '-ascii');
fprintf('       - Written tau u (v point) data to %s\n',outname);
% Tauy at u point (g_tauy_u == gowtauvu)
outname = [str(2).dir '/' str(2).exp '.tauy_u.dat'];
c = wstress(:,:,3); b = permute(c, [2 1]); a = reshape(b, [imax*jmax 1]);
save (outname, 'a', '-ascii');
fprintf('       - Written tau v (u point) data to %s\n',outname);
% Tauy at v point (g_tauy_v == gowtauvv)
outname = [str(2).dir '/' str(2).exp '.tauy_v.dat'];
c = wstress(:,:,4); b = permute(c, [2 1]); a = reshape(b, [imax*jmax 1]);
save (outname, 'a', '-ascii');
fprintf('       - Written tau v (v point) data to %s\n',outname);
%
% U wind velocity 
outname = [str(2).dir '/' str(2).exp '.wvelx.dat'];
c = wvelocity(:,:,1); b = permute(c, [2 1]); a = reshape(b, [imax*jmax 1]);
save (outname, 'a', '-ascii');
fprintf('       - Written u wind speed data to %s\n',outname);
% V wind speed
outname = [str(2).dir '/' str(2).exp '.wvely.dat'];
c = wvelocity(:,:,2); b = permute(c, [2 1]); a = reshape(b, [imax*jmax 1]);
save (outname, 'a', '-ascii');
fprintf('       - Written v wind speed data to %s\n',outname);
%
% Save 2-D ASCII wind speed scalar for BIOGEM
outname = [str(2).dir '/' str(2).exp '.windspeed_uvaa.dat'];
a = wspeed_uvaa;
save(outname,'a','-ascii');
outname = [str(2).dir '/' str(2).exp '.windspeed_uvma.dat'];
a = wspeed_uvma;
save(outname,'a','-ascii');
fprintf('       - Written BIOGEM windspeed data to %s\n',outname);
% Save 2-D ASCII wind speed scalar on FULL grid 
%   no seperation between land-ocean while regridding
%   uvaa regrids annual average wind velocities to derive speed
wspeed_all = gowspdall;
outname = [str(2).dir '/' str(2).exp '.windspeed_uvaa_all.dat'];
a = wspeed_all;
save(outname,'a','-ascii');
fprintf('       - Written FULL GRID uvaa windspeed data\n');
% Save 2-D ASCII wind speed scalar on FULL grid 
%   no seperation between land-ocean while regridding
%   uvma regrids monthly wind velocities to derive annual mean speed
wspeed_all = gowspdall_uvma;
outname = [str(2).dir '/' str(2).exp '.windspeed_uvma_all.dat'];
a = wspeed_all;
save(outname,'a','-ascii');
fprintf('       - Written FULL GRID uvma windspeed data\n');
% Save 2-D ASCII wind speed scalar on ocean + land grids
%   only consider GCM ocean grids for cGENIE ocean windspeed
%   only consider GCM land grids for cGENIE land windspeed
wspeed_ents = wspeed_uvaa+wspeed_uvaal; % combine ocean+land
outname = [str(2).dir '/' str(2).exp '.windspeed_uvaa_ents.dat'];
a = wspeed_ents;
save(outname,'a','-ascii');
fprintf('       - Written ENTS uvaa windspeed (derived from annual mean velocity)\n');
% Save 2-D ASCII wind speed scalar on ocean + land grids
%   only consider GCM ocean grids for cGENIE ocean windspeed
%   only consider GCM land grids for cGENIE land windspeed
wspeed_ents = wspeed_uvma+wspeed_uvmal; % combine ocean+land
outname = [str(2).dir '/' str(2).exp '.windspeed_uvma_ents.dat'];
a = wspeed_ents;
save(outname,'a','-ascii');
fprintf('       - Written ENTS uvma windspeed (derived from monthly velocity)\n');

%
% *********************************************************************** %

% *********************************************************************** %
% *** END *************************************************************** %
% *********************************************************************** %
%
%
% *********************************************************************** %
