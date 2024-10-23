function [grid_borders,cell_borders] = find_grid_borders_update(grid_borders,cell_borders,grid_islands,grid_mask,n_islands)
%
%%

% *********************************************************************** %
% *** REFINE ISLANDS BORDERS ******************************************** %
% *********************************************************************** %
%
% determine mask size (remember: [rows columns])
[jmax, imax] = size(grid_islands);
% create search array
vdsrch = [1 1; 1 0; 1 -1; 0 -1; -1 -1; -1 0; -1 1; 0 1]; 
vdsrch_nsew = [1 0; 0 1; -1 0; 0 -1]; 
% copy & expand grids
% NOTE: mask is defined with '1' for ocean
gb_ex = grid_borders;
gb_ex = [gb_ex(:,end) gb_ex gb_ex(:,1)];
gb_ex = [gb_ex(1,:); gb_ex; gb_ex(end,:)];
gb_ex(1,:)   = 0;
gb_ex(end,:) = 0;
gi_ex = grid_islands;
gi_ex = [gi_ex(:,end) gi_ex gi_ex(:,1)];
gi_ex = [gi_ex(1,:); gi_ex; gi_ex(end,:)];
gi_ex(1,:)   = 0;
gi_ex(end,:) = 0;
gm_ex = grid_mask;
gm_ex = [gm_ex(:,end) gm_ex gm_ex(:,1)];
gm_ex = [gm_ex(1,:); gm_ex; gm_ex(end,:)];
gm_ex(1,:)   = 0;
gm_ex(end,:) = 0;
% find any non zero island number cell bordering poles
% and set extended array accordingly
if (max(grid_islands(1,:)) > 0)
    gi_ex(1,:) = max(grid_islands(1,:));
elseif (min(grid_islands(1,:)) == -3)
    gi_ex(1,:) = -3;
end
if (max(grid_islands(end,:)) > 0)
    gi_ex(end,:) = max(grid_islands(end,:));
elseif (min(grid_islands(end,:)) == -3)
    gi_ex(1,:) = -3;
end
%
% *** ASSIGN BORDER NUMBER ********************************************** %
%
% search across grid -- assign border number
for i = 2:imax+1
    for j = 2:jmax+1
        % test for cell being a border
        if gb_ex(j,i)
            % set results index empty
            isrch = 0;
            % search surrounding cells
            for s = 1:length(vdsrch)
                loc_j = j + vdsrch(s,1);
                loc_i = i + vdsrch(s,2);
                % test for adjacent island cell
                if (gi_ex(loc_j,loc_i) ~= 0)
                    isrch = gi_ex(loc_j,loc_i);
                    % set border number in cell array if that particular
                    % island number has not been recorded yet
                    % NOTE: remember to convert back array indices
                    if (isempty(find(cell_borders{j-1,i-1} == isrch)))
                        cell_borders{j-1,i-1} = [cell_borders{j-1,i-1} isrch];
                    end
                end
            end
            % set border number -- original grid array
            % NOTE: remember to convert back array indices
            grid_borders(j-1,i-1) = isrch;
        end
    end
end
%
% *** CLEAN UP HANGING BORDER CELLS ************************************* %
%
% search across grid -- remove any border cells surrounded on 3 sides by 
% the same island and/or a pole cell
% this should prevent the double-back issue wihtout need for
% post-processing the partially created path in make_paths
% NOTE: cell_borders = cell([jmax imax]) (NOT on extended grid)
% search all islands
% NOTE: iterative search
for islnd = 1:n_islands
    search = true;
    while search
        search = false;
        for i = 2:imax+1
            for j = 2:jmax+1
                % test for cell being a border and record cell element index
                % test whether border cell has a grid point from the same island on EACH side
                loc_islndi = find(cell_borders{j-1,i-1} == islnd);
                if (~isempty(loc_islndi))
                    if ( (gi_ex(j-1,i) == islnd) && (gi_ex(j+1,i) == islnd) ) || ( (gi_ex(j,i-1) == islnd) &&  (gi_ex(j,i+1) == islnd) )
                        % mark current cell as an island for the purpose of the search
                        % (even though this is not an island)
                        gi_ex(j,i) = islnd;
                        % mark border cell element
                        cell_borders{j-1,i-1}(loc_islndi) = 0;
                        % update old format border array for displayy purposes
                        % NOTE: not extended grid
                        grid_borders(j-1,i-1) = 0;
                        % allow search to continue
                        search = true;
                        % report action
                        disp(['         -> Border #' num2str(islnd) ' grid point removed at @ (' num2str(i-1) ',' num2str(jmax-(j-1)+1) ')']);
                    end
                end
            end
        end
    end
end
%
% *********************************************************************** %
% *** END *************************************************************** %
% *********************************************************************** %
