% cookiegen_settings

% *********************************************************************** %
% *** USER SETTINGS FOR cookiegen CONFIG GENERATOR ********************** %
% *********************************************************************** %
%
% *** CONFIG NAME AND INPUT DATA SETTINGS ******************************* %
%
par_wor_name='wordrake';       % ['STRING'] 8-char (output) config name
par_gcm='mask';                % ['STRING'] input format/GCM name
par_expid='wordrake';          % ['STRING'] input experiment/data name
par_age=0.0;                   % [0.0-4570.0] optional age (Myr)
%
% *** INPUT + OUTPUT SETTINGS ******************************************* %
%
par_pathin='INPUT';            % ['STRING'] path to input dir
par_pathout='OUTPUT';          % ['STRING'] path to output dir
par_plotformat='png';          % 'png' ,'jpg','pdf' (high res), '' for NONE
opt_user=true;                 % [false/true] force user input to grid
%
% *** GRID -- HORIZONTAL ************************************************ %
%
par_max_i=36;                  % [1-72] # grid cells in longitude dir (i)
par_max_j=36;                  % [1-72] # grid cells in latitude  dir (j)
opt_equalarea=true;            % [false/true] equal area grid?
par_lon_off=-180.0;            % [-360-0] longitude offset of grid start
%
% *** GRID -- vertical ************************************************** %
%
par_max_k=16;                  % [1-99] total # levels in ocean
par_max_k_shallow=15;          % [1-99] maximum shallow water # level
par_min_k=3;                   % [1-99] minimum ocean depth (k value)
par_max_D=5000.0;              % [0.0-99999.9] ocean depth (m) at k=1
%
% *** MODULE SUPPORT **************************************************** %
%
opt_makegold=true;             % [false/true] make GOLDSTEIN ocean files?
opt_makeseds=false;            % [false/true] make SEDGEM files
opt_makeents=false;            % [false/true] make ENTS files?
%
% *** BOUNDARY CONDITION SIMPLIFICATIONS ******************************** %
%
opt_makezonalwind=true;       % [false/true] generate zonal winds
opt_makezonalalbedo=true;     % [false/true] generate zonal average albedo
opt_makerndseds=false;        % [false/true] generate randomised sed depths
%
% *********************************************************************** %
