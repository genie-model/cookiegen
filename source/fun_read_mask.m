function [mask] = fun_read_mask(str)
%
%%

% *********************************************************************** %
% *** READ AND RETURN MASK ********************************************** %
% *********************************************************************** %
%
% *** READ IN MASK FILE ************************************************* %
%
% Create filename
if isempty(str(1).path)
    loc_str_file = [str(1).mask '.txt'];
else
    loc_str_file = [str(1).path '/' str(1).mask '.txt'];
end
% Check for file existing ... :o)
% ... then load.
if (exist(loc_str_file, 'file') == 2)
    loc_mask = load(loc_str_file);
else
    disp(['       ERROR: mask file: ', str(1).mask, ' does not exist (at location ' str(1).path ').']);
    disp(['              Problem with correct extension? (mask file extension must be .txt)']);
    disp(['--------------------------------------------------------']);
    disp([' ']);
    return;
end
%
% *** DERIVE MASK ******************************************************* %
%
% (for compatability only)
gm = loc_mask;
%
% *** RETURN ARRAYS ***************************************************** %
%
mask = gm;
%
% *********************************************************************** %
% *** END *************************************************************** %
% *********************************************************************** %
