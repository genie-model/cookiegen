function [] = fun_eval_bath(PLONI,PLATI,PZI,PLONO,PLATO,PZO)
% fun_eval_bath
%
%   ***********************************************************************
%   *** evaluate bathymetry profile ***************************************
%   ***********************************************************************
%
%
%   ***********************************************************************
%   *** HISTORY ***********************************************************
%   ***********************************************************************
%
%   26/09/07: CREATED
%
%   ***********************************************************************

% *********************************************************************** %
% *** INITIALIZE PARAMETERS & VARIABLES ********************************* %
% *********************************************************************** %
%
% process dummy parameters
loni = PLONI;
lati = PLATI;
lono = PLONO;
lato = PLATO;
gzi  = PZI;
gzo  = PZO;
%
% *** misc (local) parameters ******************************************* %
%
% constants
par_rEarth = 6371000.0;
par_yrtos = 365.25*24.0*3600.0;
% set date
str_date = [datestr(date,11), datestr(date,5), datestr(date,7)];
%
% *** create grids and dimensions *************************************** %
%
% create 2D lat & lon grids
[gi_latn gi_lone] = meshgrid(lati(2:end),loni(2:end));
[gi_lats gi_lonw] = meshgrid(lati(1:end-1),loni(1:end-1));
[go_latn go_lone] = meshgrid(lato(2:end),lono(2:end));
[go_lats go_lonw] = meshgrid(lato(1:end-1),lono(1:end-1));
% area
gi_area(:,:) = 2.0*pi*(par_rEarth^2)*(sin((pi/180.0)*gi_latn) - sin((pi/180.0)*gi_lats)).*((gi_lone-gi_lonw)/360.0);
go_area(:,:) = 2.0*pi*(par_rEarth^2)*(sin((pi/180.0)*go_latn) - sin((pi/180.0)*go_lats)).*((go_lone-go_lonw)/360.0);
%
% *********************************************************************** %

% *********************************************************************** %
% *** EVAL BATH ********************************************************* %
% *********************************************************************** %
%
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

% here, we could re-weight grid area 
% depending on the topographic height of the grid cell
% note that gi_z remains as height and not depth

% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%
% *** create hypsographic curve for input grid ************************** %
% NOTE: go depths regridded from only abyssal gi cells
%       are not exactly the same value ... allow for 1% uncertainty
%
% find the different depth layers in the original grid
vzi_nonan = gzi(find(~isnan(gzi)));
zi = unique(vzi_nonan);
% calculate total seafloor area occupied by each depth
for z=1:length(zi)
    zi_area(z) = sum(gi_area(find(gzi==zi(z))),'all');
end
% calculate total seafloor area occupied by each zi depth
for z=length(zi):-1:1
    zi_cumarea(z) = sum(gi_area(find(gzi>=1.01*zi(z))),'all');
end
% 
% *** create hypsographic curve for output grid ************************* %
%
% % calculate total seafloor area occupied by each depth
% for z=1:length(zi)
%     zo_area(z) = sum(go_area(find(gzo==zi(z))),'all');
% end
% calculate total seafloor area occupied by each zi depth
for z=length(zi):-1:1
    zo_cumarea(z) = sum(go_area(find(gzo>=1.01*zi(z))),'all');
end
%
% *** plot hypsographic curves ****************************************** %
%
hold on;
scatter(zi_area/1.0E12,zi,50.0,'r','filled');
plot(zi_cumarea/1.0E12,zi,'r','LineWidth',1.0);
% scatter(zo_area/1.0E12,zi,50.0,'b','filled');
plot(zo_cumarea/1.0E12,zi,'b--','LineWidth',2.0);
axis([0 4.5E14/1.0E12 -6000 0]);
xlabel('Area (million km^2)');
ylabel('Topographic height (m)');
legend({'Area of each input depth level','Cumulative area (input)','Cumulative area (regridded)'})
%
% *********************************************************************** %

% *********************************************************************** %
% *** END *************************************************************** %
% *********************************************************************** %
%
% return hypsographic curve data?
%
% *********************************************************************** %
