function [n_paths,v_paths,n_islands,grid_paths] = make_paths(cell_borders,n_islands,i_poles);
%
%%

% *********************************************************************** %
% *** ELUCIDATE PATHS *************************************************** %
% *********************************************************************** %
%
% NOTE: polar borders can be incomplete, hence the pole must be followed
%       in the same way as an island is
%       (special borders are created to enable this to occur automatically)
% NOTE: don't worry initially about the rows being counted from
%       top-top-bottom (opposite of GENIE grid)
%       => this is fixed at the very end
% NOTE: format is of first parameter being the direction to the *next*
%       cell
% NOTE: path cannot double-back on itself through the same cell
% NOTE: GOLDSTEIN directions:
%        2 == North (+ve v)
%       -2 == South (-ve v)
%        1 == East  (+ve u)
%       -1 == West  (-ve u)
% determine mask size (remember: [rows columns])
[jmax imax] = size(cell_borders);
% create search array + extended array
% NOTE: remember (j,i) and rows are flipped (South is up)
% directions: 1==S, 2==E, 3==N, 4==W
vdsrch    = [1 0; 0 1; -1 0; 0 -1];
vdsrch_ex = [vdsrch; vdsrch];
% expand cell array into n-dimensonal grid borders array
% NOTE: assume additional dimension will be the same as n_islands
grid_borders_n = zeros(jmax,imax,n_islands);
for j = 1:jmax
    for i = 1:imax
        for islnd = 1:n_islands
            loc_k = find(cell_borders{j,i}==islnd);
            if ~isempty(loc_k)
                grid_borders_n(j,i,islnd) = islnd;
            end
        end
    end
end
% copy & expand border grid
gbn_ex = grid_borders_n;
gbn_ex = [gbn_ex(:,end,:) gbn_ex gbn_ex(:,1,:)];
gbn_ex = [gbn_ex(1,:,:); gbn_ex; gbn_ex(end,:,:)];
gbn_ex(1,:,:)   = 0.0;
gbn_ex(end,:,:) = 0.0;
% created 'searched' grid
gsn_ex = zeros(jmax+2,imax+2,n_islands);
% initialize paths arrays
v_paths = [];
n_paths = [];
% create naked polar island borders
% NOTE: count how many polar borders have been created ...
n_poles = 0;
if (~isempty(i_poles))
    % N pole
    if ~isempty(find(i_poles == -1))
        n_islands = n_islands + 1;
        gbn_ex(2,:,n_islands) = n_islands;
        n_poles = n_poles + 1;
    end
    % S pole
    if ~isempty(find(i_poles == -2))
        n_islands = n_islands + 1;
        gbn_ex(end-1,:,n_islands) = n_islands;
        n_poles = n_poles + 1;
    end
end
%
if (n_islands >= 2)
    %
    disp(['       * Ignoring border #1']);
    % set dummy first path data -- length 1
    v_paths = [v_paths; 0 0 0];
    n_paths = [n_paths 1];
    % LOOP >>>
    for islnd = 2:n_islands
        %
        disp(['       * Creating raw path #' num2str(islnd) ' ...']);
        % find upper LH border cell of the current island border
        % => must be a uppermost border cell with ocean to the West
        %    EXCEPT for a polar island 
        %    (ignoring what happens if land stretches 360 along Eq.)
        % NOTE: for an island connecte to the N pole, the corner found
        %       will not be the uppermost corner
        % raster left-to-right, top-to-bottom
        corner = false;
        for j = 2:jmax+1
            for i = 2:imax+1
                if (islnd <= (n_islands - n_poles))
                    if ((gbn_ex(j,i,islnd) == islnd) && (gbn_ex(j,i-1,islnd) == 0)), corner = true; end
                else
                    if (gbn_ex(j,i,islnd) == islnd), corner = true; end
                end
                if corner, break; end
            end
            if corner, break; end
        end
        % check upper LH border cell was even found!
        if (~corner)
            diary off;
            error(['Error. \nCould not find upper LH corner of border #' num2str(islnd) ' : %s'],'Exiting ...');
            return;
        end
        % copy initial location
        loc_j = j;
        loc_i = i;
        % mark initial cell as searched
        gsn_ex(loc_j,loc_i,islnd) = 1;
        % test for E-W wall and also mark wrap-around cell as searched
        if (loc_i == 2), gsn_ex(loc_j,imax+2,islnd) = 1; end
        if (loc_i == imax+1), gsn_ex(loc_j,1,islnd) = 1; end
        % check assumption of moving East being a valid direction
        % NOTE: in reporting location, remember that the grid is extended 
        %       and transposed (but leave up-side-down)
        if (gbn_ex(loc_j,loc_i+1,islnd) ~= islnd)
            diary off;
            error(['Error. \nFailed to follow border #' num2str(islnd) ' in E direction starting from: (',num2str(loc_j-1),',',num2str(loc_i-1),'): %s'],'Exiting ...');
            return;
        end
        % initial direction to the next cell East [1]
        % => take first step in that direction
        % record direction and current location (in core grid indices)
        % NOTE: 1 == East in GOLDSTEIN path notation
        v_paths = [v_paths; 1 loc_i-1 loc_j-1];
        % now move 1 East
        loc_i = loc_i + 1;
        % test for E-W wall ...
        % mark i==1 cell as implicitly, already searched
        if (loc_i == imax+2)
            loc_i = 2;
            gsn_ex(loc_j,loc_i,islnd)  = 1;
            gsn_ex(loc_j,imax+2,islnd) = 1;
        end
        % mark as searched
        gsn_ex(loc_j,loc_i,islnd) = 1;
        % record direction taken
        loc_s = 2;
        % initialize vector length at 1
        % (as the vector has already been populated with its first line)
        n_path = 1;
        % now follow path around island -- clockwise
        % NOTE: for an island connected to the N pole, the direction
        %       will be ANTICLOCKWISE (island on left) 
        follow = true;
        while follow
            % search surrounding cells ... 
            % ... find adajacent, unmarked border cell
            % => if no unmarked border cells exist, finish ...
            %    ... but ONLY if the marked border cell is the start
            %    otherwise, the path has doubled back on itself
            %    and a segment of the path will need to be removed
            % NOTE: search direction should start to the RIGHT of last move
            %       direction and change search direction anticlockwise
            %       (keeping the island to the right at all times)
            follow = false;
            for s = loc_s-1:loc_s+2
                loc_jj = loc_j + vdsrch_ex(s,1);
                loc_ii = loc_i + vdsrch_ex(s,2);
                % test for adjacent border cell that is
                % the current border AND unmarked
                if (gbn_ex(loc_jj,loc_ii,islnd) == islnd) && ~gsn_ex(loc_jj,loc_ii,islnd)
                    % record:
                    % *** current location ***
                    % AND 
                    % *** direction to next cell ***
                    % directions: 1==S, 2==E, 3==N, 4==W
                    % GOLDSTEIN path notation reminder:
                    %        2 == North
                    %       -2 == South
                    %        1 == East
                    %       -1 == West
                    % NOTE: take into account expanded i,j indices and
                    %       record location in core grid coordinates
                    switch s
                        case {1,5}
                            v_paths = [v_paths; -2 loc_i-1 loc_j-1];
                        case {2,6}
                            v_paths = [v_paths;  1 loc_i-1 loc_j-1];
                        case {3,7}
                            v_paths = [v_paths;  2 loc_i-1 loc_j-1];
                        case {4,8}
                            v_paths = [v_paths; -1 loc_i-1 loc_j-1];
                    end
                    % update path length count
                    n_path = n_path + 1;
                    % copy location
                    loc_j = loc_jj;
                    loc_i = loc_ii;
                    % mark tested cell as searched
                    gsn_ex(loc_j,loc_i,islnd) = 1;
                    % test for E-W wall:
                    % adjust (j,i) location if necessary and mark as searched
                    if (loc_ii == 1)
                        gsn_ex(loc_j,imax+1,islnd) = 1; % wrap-around location
                        loc_i = imax+1;
                        gsn_ex(loc_j,imax+2,islnd) = 1; % implicitly, already searched
                    elseif (loc_ii == 2)
                        gsn_ex(loc_j,imax+2,islnd) = 1; % wrap-around location
                    elseif (loc_ii == imax+2)
                        gsn_ex(loc_j,2,islnd) = 1;      % wrap-around location
                        loc_i = 2;
                        gsn_ex(loc_j,1,islnd) = 1;      % implicitly, already searched
                    elseif (loc_ii == imax+1)
                        gsn_ex(loc_jj,1,islnd) = 1;     % wrap-around location
                    end
                    % record direction taken
                    % NOTE: remember that s might be > 4 and in the
                    %       wrap-around part of the search vector
                    %       or 1, when loc_s-1 in the next search loop
                    %       will casue a problem ...
                    if (s > 4)
                        loc_s = s - 4;
                    else
                        loc_s = s;
                    end
                    if (loc_s == 1), loc_s = loc_s + 4; end
                    % continue ...
                    follow = true;
                    % exit (s) loop
                    break;
                end % end if
            end % end s loop

            if ~follow
                % at this point, no progress on un-searched path is possible
                % => find first search cell -- this should be the START ...
                % ... ...
                %
                % find FIRST border cell that is current border AND unmarked
                % NOTE: remember (i,j) is START location on the extended grid
                % NOTE: loc_s is not adjusted when the path search draws a blank
                for s = loc_s-1:loc_s+2
                    loc_jj = loc_j + vdsrch_ex(s,1);
                    loc_ii = loc_i + vdsrch_ex(s,2);
                    if (gbn_ex(loc_jj,loc_ii,islnd) == islnd) && gsn_ex(loc_jj,loc_ii,islnd)
                        if ((loc_jj == j) && ((loc_ii == i) || (loc_ii-imax == i)))
                            % start cell found -- true end of path!
                            % exit (s) loop
                            break;
                        else
                            % doubling back situation ... :(
                            % find last time this cell was serched
                            % NOTE: n_path is the local (current) path legnth
                            %       v_paths contains all path elements + 1
                            loc_v_path = [];
                            [loc_nmax tmp] = size(v_paths);
                            for n=[loc_nmax:-1:loc_nmax-n_path]
                                if (v_paths(n,2) == loc_ii-1) && (v_paths(n,3) == loc_jj-1)
                                    loc_v_path = v_paths(n,:);
                                    break; % exit (n) loop
                                end
                            end
                            %
                            if isempty(loc_v_path)
                                disp([' *** There is a very very very sad issue somehow :( ']);
                                disp([' ']);
                                diary off;
                                return;
                            end
                        end
                        % move (return) to duplicate cell location
                        % and test for E-W wall
                        loc_i = v_paths(n,2) + 1;
                        loc_j = v_paths(n,3) + 1;
                        if (loc_ii == 1)
                            loc_i = imax+1;
                        elseif (loc_ii == imax+2)
                            loc_i = 2;
                        end
                        % remove path back to and including dup cell
                        v_paths(n:loc_nmax,:) = [];
                        % update path length
                        n_path = n_path - (loc_nmax-n+1);
                        % set search direction (East of the last move direction)
                        % directions: 1==S, 2==E, 3==N, 4==W
                        % GOLDSTEIN path notation reminder:
                        %        2 == North
                        %       -2 == South
                        %        1 == East
                        %       -1 == West
                        switch v_paths(end,1)
                            case {2}
                                loc_s = 3;
                            case {-2}
                                loc_s = 5;
                            case {1}
                                loc_s = 2;
                            case {-1}
                                loc_s = 4;
                        end
                        % continue ...
                        follow = true;
                        % exit (s) loop
                        break;
                    end
                end
            end

        end % end follow while




        % add final location, calculate direction to start, update count
        % NOTE: take into account upside-down GENIE array MATLAB indexing
        % reminder:
        %        2 == North
        %       -2 == South
        %        1 == East
        %       -1 == West
        if ((loc_i-i+1) == imax) || ((loc_i-i) == -1) %East
            v_paths = [v_paths;  1 loc_i-1 loc_j-1];
        elseif ((loc_i-i-1) == -imax) || ((loc_i-i) == 1) %West
            v_paths = [v_paths; -1 loc_i-1 loc_j-1];
        elseif (loc_j-j) == 1 %North
            v_paths = [v_paths;  2 loc_i-1 loc_j-1];
        elseif (loc_j-j) == -1 %South
            v_paths = [v_paths; -2 loc_i-1 loc_j-1];
        else
            disp([' *** There is a sad issue somehow :( ']);
            disp([' ']);
            diary off;
            error(['Error. \nFailed to complete path loop @ (',num2str(loc_i),',',num2str(jmax-loc_j+1),'): %s'],'Exiting ...');
            return;
        end
        %
        n_path = n_path + 1;
        % write out path length
        n_paths = [n_paths n_path];
        %
    end
    % <<< LOOP
end
% now fix the up-side-down World
for p = 1:sum(n_paths)
    v_paths(p,3) = jmax-v_paths(p,3)+1;
end
%
% *** PATH POST-PROCESSING ********************************************** %
%
% Included in the automatically generated paths, are 2 illegal moves that
% result in the c-grid being transversed twice
% (and hence the path direction ambigeous).
% There are:
% The top right corner cell of a (bend in a) path,
% whether an outer or innter bend, and characterized by:
% path cells to the E and the S
% (and ocean and/or land to the W and the N)
% Both other corner configurations are valid.
%
% The search is basically for:
% (a) a '+1' followed by a '-2'
% (b) a '+2' followed by a '-1'
% (and remembering to deal with the wrap-around E-W grid)
%
% NOTE: don't bother replacing the 1st path (not written out!)
%
% REMEMBER:
%
%   ----v----
%   |       |
%   |   c   u
%   |       |
%   ---------
%
%        2 == North (+ve v)
%       -2 == South (-ve v)
%        1 == East  (+ve u)
%       -1 == West  (-ve u)
%
if n_islands >= 2
    %
    for islnd = 2:n_islands
        %
        disp(['       * Building path #' num2str(islnd) ' ...']);
        %
        pmin = sum(n_paths(1:islnd-1))+1;
        pmax = sum(n_paths(1:islnd));
        dp = pmax-pmin+1;
        % create tempoary wrap-around array
        % (in terms of following the path loop)
        vp_ex = v_paths(pmin:pmax,:);
        vp_ex = [vp_ex(end,:); vp_ex; vp_ex(1,:)];
        %
        vpp = [];
        %
        p = 2;
        pp = 1;
        while p <= dp+1
            if (vp_ex(p-1,1)==1 && vp_ex(p,1)==-2) || (vp_ex(p-1,1)==2 && vp_ex(p,1)==-1)
                % top RH corner, clockwise
                % [p-1] == 1; [p] == -2
                % BORDER DIRECTIONS   GOLDSTEIN VELOCITY
                % -----------------   ----------v------
                % !-2 >-1 > p ! o !   !-2 >-1 >   ! o !
                % ----------v------   ----------v------
                % ! X ! X ! 1 ! o !   ! X ! X ! 1 ! o !
                % ----------v------   ----------v------
                % ! X ! X ! 2 ! o !   ! X ! X ! 2 ! o !
                % ----------v------   -----------------
                % top RH corner, anticlockwise
                % [p-1] == 2; [p] == -1
                % BORDER DIRECTIONS   GOLDSTEIN VELOCITY
                % -----------------   -----------------
                % < 2 < 1 < p ! o !   ! 2 < 1 <   < o !
                % ----------^------   ----------^------
                % ! X ! X !-1 ! o !   ! X ! X !-1 ! o !
                % ----------^------   ----------^------
                % ! X ! X !-2 ! o !   ! X ! X !-2 ! o !
                % -----------------   -----------------
                disp(['         -> NE corner :: ' 'Skip path entry @ (' num2str(vp_ex(p,2)) ',' num2str(vp_ex(p,3)) ')']);
                % NOTE: add no path component (or pp count update)
                % update p count
                p = p+1;
            elseif (vp_ex(p-1,1)==-2 && vp_ex(p,1)==1) || (vp_ex(p-1,1)==-1 && vp_ex(p,1)==2)
                % bottom LH corner, anticlockwise
                % [p-1] == -2; [p] == 1
                % BORDER DIRECTIONS   GOLDSTEIN VELOCITY
                % -----------------   ------v----------
                % ! o !-2 ! X ! X !   ! o !-2 ! X ! X !
                % ------v----------   ------v----------
                % ! o !-2 ! X ! X !   ! o !-1 ! X ! X !
                % ------v----------   ------v----------
                % ! o ! p > 1 > 2 >   ! o !pp > 1 > 2 >
                % -----------------   -----------------
                % bottom LH corner, clockwise
                % [p-1] == -1; [p] == 2
                % BORDER DIRECTIONS   GOLDSTEIN VELOCITY
                % ------^----------   ------^----------
                % ! o ! 2 ! X ! X !   ! o ! 2 ! X ! X !
                % ------^----------   ------^----------
                % ! o ! 1 ! X ! X !   ! o ! 1 ! X ! X !
                % ------^----------   ------^----------
                % ! o ! p <-1 <-2 !   ! o !pp <-1 <-2 <
                % -----------------   -----------------
                disp(['         -> SE corner :: ' 'Add additional path entry @ (' num2str(vp_ex(p,2)) ',' num2str(vp_ex(p,3)) ')']);
                % duplicate path component
                % NOTE: adjust vector in first entry to complete path turn
                vpp(pp,:) = vp_ex(p,:);
                vpp(pp,1) = vp_ex(p-1,1);
                pp = pp+1;
                vpp(pp,:) = vp_ex(p,:);
                pp = pp+1;
                % update p count
                p = p+1;
            elseif (vp_ex(p-1,1)==2 && vp_ex(p,1)==1) || (vp_ex(p-1,1)==1 && vp_ex(p,1)==2) || (vp_ex(p-1,1)==1 && vp_ex(p,1)==1) || (vp_ex(p-1,1)==2 && vp_ex(p,1)==2)
                % top LH corner, clockwise
                % [p-1] == 2; [p] == 1
                % BORDER DIRECTIONS   GOLDSTEIN VELOCITY
                % -----------------   -----------------
                % ! o ! p > 1 > 2 >   ! o ! p > 1 > 2 >
                % ------^----------   ------^----------
                % ! o !-1 ! X ! X !   ! o !-1 ! X ! X !
                % ------^----------   ------^----------
                % ! o !-2 ! X ! X !   ! o !-2 ! X ! X !
                % -----------------   -----------------
                % bottom RH corner, anticlockwise
                % [p-1] == 2; [p] == 1
                % BORDER DIRECTIONS   GOLDSTEIN VELOCITY
                % ----------^------   ----------^------
                % ! X ! X !-2 ! o !   ! X ! X !-2 ! o !
                % ----------^------   ----------^------
                % ! X ! X !-1 ! o !   ! X ! X !-1 ! o !
                % ----------^------   ----------^------
                % ! 2 > 1 > p ! o !   ! 2 > 1 > p ! o !
                % -----------------   -----------------
                % East (right), North (up)
                % BORDER DIRECTIONS   GOLDSTEIN VELOCITY
                % ------^----------   ------^----------
                % ! o ! 1 ! X ! X !   ! o ! 1 ! o ! o !
                % ------^----------   ------^----------
                % !-1 > p > 1 > 1 >   !-1 > p > 1 > 2 >
                % ------^----------   ------^----------
                % ! o !-1 ! o ! o !   ! o !-1 ! o ! o !
                % -----------------   -----------------
                vpp(pp,:) = vp_ex(p,:);
                pp = pp+1;
                p = p+1;
            elseif (vp_ex(p-1,1)==-1 && vp_ex(p,1)==-2) || (vp_ex(p-1,1)==-2 && vp_ex(p,1)==-1) || (vp_ex(p-1,1)==-1 && vp_ex(p,1)==-1) || (vp_ex(p-1,1)==-2 && vp_ex(p,1)==-2)
                % top LH corner, anticlockwise
                % [p-1] == -1; [p] == -2
                % BORDER DIRECTIONS   GOLDSTEIN VELOCITY
                % -----------------   -----------------
                % ! o ! p <-1 <-2 !   ! o ! p <-1 <-2 <
                % ------v----------   ------v----------
                % ! o ! 1 ! X ! X !   ! o ! 1 ! X ! X !
                % ------v----------   ------v----------
                % ! o ! 2 ! X ! X !   ! o ! 2 ! X ! X !
                % ------v----------   -----------------
                % bottom RH corner, clockwise
                % [p-1] == -2; [p] == -1
                % BORDER DIRECTIONS   GOLDSTEIN VELOCITY
                % -----------------   ----------v------
                % ! X ! X ! 2 ! o !   ! X ! X ! 2 ! o !
                % ----------v------   ----------v------
                % ! X ! X ! 1 ! o !   ! X ! X ! 1 ! o !
                % ----------v------   ----------v------
                % <-2 <-1 < p ! o !   !-2 <-1 < p ! o !
                % -----------------   -----------------
                % West (left), South (down)
                % BORDER DIRECTIONS   GOLDSTEIN VELOCITY
                % -----------------   ------v----------
                % ! o !-1 ! X ! X !   ! o !-1 ! o ! o !
                % ------V----------   ------v----------
                % < 1 < p <-1 <-1 !   ! 1 < p <-1 <-2 <
                % ------V----------   ------v----------
                % ! o ! 1 ! o ! o !   ! o ! 1 ! o ! o !
                % ------V----------   -----------------
                vpp(pp,:) = vp_ex(p,:);
                vpp(pp,1) = vp_ex(p-1,1);
                pp = pp+1;
                p = p+1;
            else
                disp(' *** IMPOSSIBLE!');
                return;
            end
            %
        end
        % add new paths back in array
        v_paths(pmin:pmax,:) = vpp(:,:);
    end
    %
end
%
%
% *** CREATE 2D PATHS MAP *********************************************** %
%
% create empty array
grid_paths = zeros(jmax,imax);
% populate array
% NOTE: don't bother outputting 1st path
if n_islands >= 2
    for p = (n_paths(1)+1):sum(n_paths)
        grid_paths(v_paths(p,3),v_paths(p,2)) = v_paths(p,1);
    end
end
% re-orientate
grid_paths = flipud(grid_paths);
%
% *********************************************************************** %

% *********************************************************************** %
% *** END *************************************************************** %
% ***********************************************************************
%
%%%
%
% *********************************************************************** %
