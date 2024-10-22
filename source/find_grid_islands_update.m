function [grid_islands,n_islands,i_islands,i_poles]  = find_grid_islands_update(grid_islands,n_islands,i_islands,grid_borders)
%
%%

% *********************************************************************** %
% *** UPDATE THE ISLANDS COUNT ****************************************** %
% *********************************************************************** %
%
% *** INITIALIZE ******************************************************** %
%
% determine mask size (remember: [rows columns])
[jmax, imax] = size(grid_islands);
% create search array
vdsrch = [1 1; 1 0; 1 -1; 0 -1; -1 -1; -1 0; -1 1; 0 1]; 
% copy & expand grids
gb_ex = grid_borders;
gb_ex = [gb_ex(:,end) gb_ex gb_ex(:,1)];
gb_ex = [gb_ex(1,:); gb_ex; gb_ex(end,:)];
gb_ex(1,:)   = 0.0;
gb_ex(end,:) = 0.0;
gi_ex = grid_islands;
gi_ex = [gi_ex(:,end) gi_ex gi_ex(:,1)];
gi_ex = [gi_ex(1,:); gi_ex; gi_ex(end,:)];
gi_ex(1,:)   = 0.0;
gi_ex(end,:) = 0.0;
%
% *** INDENTIFY ISLANDS CONNECTED TO POLES ****************************** %
%
% search N polar row
for i = 1:imax
    % test for cell being a N island cell
    if (grid_islands(1,i) > 0)
        % update island numbers list
        % (do this before replacing the critical number in question  ...)
        i_islands(find(i_islands == grid_islands(1,i))) = [];
        % update islands count
        n_islands = n_islands - 1;
        % replace island number
        grid_islands(find(grid_islands == grid_islands(1,i))) = -1;
    end
end
% search S polar row
for i = 1:imax
    % test for cell being a S island cell
    if (grid_islands(jmax,i) > 0)
        % update island numbers list
        % (do this before replacing the critical number in question  ...)
        i_islands(find(i_islands == grid_islands(jmax,i))) = [];
        % update islands count
        n_islands = n_islands - 1;
        % replace island number
        grid_islands(find(grid_islands == grid_islands(jmax,i))) = -2;
    elseif (grid_islands(jmax,i) == -1)
        % create a special index ...
        grid_islands(find(grid_islands == -1)) = -3;
    end
end
%
% *** RENUMBER ********************************************************** %
%
% RENUMBER NON-POLAR ISLANDS, THEN ADD BACK POLAR CONNECTED ISLANDS
% NOTE: at this point, polar-only islands are not identified
i_poles = [];
% renumber islands -- non-polar
if (n_islands > 0)
    for islnd = 1:n_islands,
        grid_islands(find(grid_islands == i_islands(islnd))) = islnd;
        i_islands(islnd) = islnd;
        disp(['       * updated find: ' num2str(islnd) ' islands']);
    end
end
% deal with poles ...
if isempty(find(grid_islands == -3))
    % renumber islands -- polar -- N
    if ~isempty(find(grid_islands == -1))
        n_islands = n_islands + 1;
        grid_islands(find(grid_islands == -1)) = n_islands;
        disp(['       * updated find: ' num2str(n_islands) ' islands']);
    else
        i_poles = [i_poles -1];
        disp(['       * updated find: N pole island (currently uncounted)']);
    end
    % renumber islands -- polar -- S
    if ~isempty(find(grid_islands == -2))
        n_islands = n_islands + 1;
        grid_islands(find(grid_islands == -2)) = n_islands;
        disp(['       * updated find: ' num2str(n_islands) ' islands']);
    else
        i_poles = [i_poles -2];
        disp(['       * updated find: S pole island (currently uncounted)']);
    end
else
    % include in the count, a double polar-connected island
    n_islands = n_islands + 1;
    grid_islands(find(grid_islands == -3)) = n_islands;
    i_poles = [i_poles -3];
    disp(['       * updated find: N-S connected polar island']);
end
%
disp(['       * total # true islands = ' num2str(n_islands-1 + length(i_poles))]);
%
% *** REORDER *********************************************************** %
%
% reorder islands such that the very largest one (i.e. most cells)
% appears first
% (ideally, this would be the one with the longest path, but may not be)
% NOTE: search (and swap) only if there is more than one island
%       (by default, a single 'real' island will be numbered first before
%       the polar island, and hence its path not included in .paths)
v_size = [];
if (n_islands > 1)
    for islnd = 1:n_islands
        n_size = length(find(grid_islands == islnd));
        v_size = [v_size n_size];
    end
    % re-order if the largest is not already #1
    % (use -3 as a temp index for swapping over indices)
    % NOTE: deal with case of >1 largest island:
    % %     take the first element always -- nmax_size = nmax_size(1);
    nmax_size = find(v_size == max(v_size));
    nmax_size = nmax_size(1); 
    if (nmax_size ~= 1)
        grid_islands(find(grid_islands == 1))         = -3;
        grid_islands(find(grid_islands == nmax_size)) = 1;
        grid_islands(find(grid_islands == -3))        = nmax_size;
    end
end
%
% *********************************************************************** %
% *** END *************************************************************** %
% *********************************************************************** %
