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
% *** CLEAN UP BURIED BORDERS ******************************************* %
%
% % NOTE: 'buried' here means
% %       (a) a border cell sandwitched between land and same # border cells
% %           and not adjacent to open (non same border) ocean
% %       (b) a border cell entirely surrounded by same # border cells
% % NOTE: treat poles as OCEAN for this (i.e. different from other filtering)
% %
% % update extended borders array
% gb_ex = grid_borders;
% gb_ex = [gb_ex(:,end) gb_ex gb_ex(:,1)];
% gb_ex = [gb_ex(1,:); gb_ex; gb_ex(end,:)];
% gb_ex(1,:)   = 0.0;
% gb_ex(end,:) = 0.0;
% % adjust ocean mask poles
% gm_ex(1,:)   = 1;
% gm_ex(end,:) = 1;
% % search across grid -- remove buried border cells
% for i = 2:imax+1
%     for j = 2:jmax+1
%         % test for cell being a border
%         if gb_ex(j,i)
%             % set island index
%             iislnd = gb_ex(j,i);
%             isrch  = 0;
%             % search surrounding cells
%             for s = 1:length(vdsrch)
%                 loc_j = j + vdsrch(s,1);
%                 loc_i = i + vdsrch(s,2);
%                 % test adjacent cells
%                 % (if ocean and not the same island border)
%                 if (gm_ex(loc_j,loc_i) && (gb_ex(loc_j,loc_i) ~= iislnd)),
%                     isrch = 1;
%                 end
%             end
%             % reset border if no adjoining ocean grid point
%             % of different border or no border exists
%             % NOTE: remember to convert back array indices
%             if ~isrch,
%                 grid_borders(j-1,i-1) = 0;
%             end
%         end
%     end
% end
%
% *** CLEAN UP LOOPS **************************************************** %
%
% -------------------------
% !   !   !   ! X !   !   ! ...
% -------------------------
% !   ! X ! X ! X ! X ! X ! ...
% -------------------------
% !   ! X !   ! X !   !   ! ...
% -------------------------
% !   ! X ! X ! X !   !   ! ...
% -------------------------
%
% -------------------------
% !   !   !   ! X !   !   ! ...
% -------------------------
% !   ! X ! X ! X !   !   ! ...
% -------------------------
% !   ! X !   ! X ! X ! X ! ...
% -------------------------
% !   ! X ! X ! X !   !   ! ...
% -------------------------
%
% % update extended borders array
% gb_ex = grid_borders;
% gb_ex = [gb_ex(:,end) gb_ex gb_ex(:,1)];
% gb_ex = [gb_ex(1,:); gb_ex; gb_ex(end,:)];
% gb_ex(1,:)   = 0.0;
% gb_ex(end,:) = 0.0;
% % revert ocean mask poles
% gm_ex(1,:)   = 1;
% gm_ex(end,:) = 1;
% % determine number of N-S-E-W connecting border cells
% gn = find_grid_ncells(grid_borders,0);
% %
% for i = 2:imax+1
%     for j = 2:jmax+1
%         % test for cell being ocean and not border
%         if (gm_ex(j,i) && (gb_ex(j,i) == 0)),
%             isrch  = 0;
%             visrch = [];
%             % search surrounding cells
%             for s = 1:length(vdsrch)
%                 loc_j = j + vdsrch(s,1);
%                 loc_i = i + vdsrch(s,2);
%                 % test adjacent cells
%                 if gb_ex(loc_j,loc_i) > 0, 
%                     isrch  = isrch + 1; 
%                     visrch = [visrch gb_ex(loc_j,loc_i)];
%                 end
%             end
%             % remove border cells with only 2 connections from loop
%             % NOTE: only if the SAME border #
%             % NOTE: catch lon grid edges
%             if isrch == 8,
%                 if (min(visrch) == max(visrch));
%                     for s = 1:length(vdsrch)
%                         loc_j = j + vdsrch(s,1);
%                         loc_i = i + vdsrch(s,2);
%                         if (loc_i == imax+2), loc_i = 2; end
%                         if (loc_i == 1), loc_i = imax+1; end
%                         if gn(loc_j-1,loc_i-1) <= 2,
%                             grid_borders(loc_j-1,loc_i-1) = 0;
%                         end
%                     end
%                 end
%             end
%         end
%     end
% end
%
% *** CLEAN UP ORPHAN BORDER CELLS ************************************** %
%
% % search across grid -- remove any orphan cells
% % NOTE: iterative search
% % initialize search as to take place
% search = true;
% while search
%     %
%     % update extended borders array
%     gb_ex = grid_borders;
%     gb_ex = [gb_ex(:,end) gb_ex gb_ex(:,1)];
%     gb_ex = [gb_ex(1,:); gb_ex; gb_ex(end,:)];
%     gb_ex(1,:)   = 0.0;
%     gb_ex(end,:) = 0.0;
%     % set search to false so only a single search ... for now
%     search = false;
%     for i = 2:imax+1
%         for j = 2:jmax+1
%             % test for cell being a border
%             if gb_ex(j,i) > 0
%                 % test for 1 or fewer N-S-E-W adjoinging border cells
%                 isrch  = 0;
%                 % search surrounding cells
%                 for s = 1:length(vdsrch_nsew)
%                     loc_j = j + vdsrch_nsew(s,1);
%                     loc_i = i + vdsrch_nsew(s,2);
%                     % test adjacent cells
%                     if gb_ex(loc_j,loc_i) > 0, isrch = isrch + 1; end
%                 end
%                 if isrch <= 1
%                     grid_borders(j-1,i-1) = 0;
%                     % found something, so try the entire search once more
%                     search = true;
%                 end
%             end
%         end
%     end
%     %
% end
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
                        disp(['         -> Element of border #' num2str(islnd) ' removed at @ (' num2str(i-1) ',' num2str(jmax-(j-1)+1) ')']);
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
