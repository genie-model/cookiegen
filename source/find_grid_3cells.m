function [grid_cells]  = find_grid_3cells(grid_mask)
%
%%

% *********************************************************************** %
% *** FIND 3 SURROUNDING LAND* CELLS ************************************ %
% *********************************************************************** %
%
%
% REMEMBER: 1 == ocean
% determine mask size (remember: [rows columns])
[jmax imax] = size(grid_mask);
gc = zeros(jmax,imax);
% add boundaries to mask
% NOTE: assign poles as land (== 1)
gm = grid_mask;
gm = [gm(:,end) gm gm(:,1)];
gm = [gm(1,:); gm; gm(end,:)];
gm(1,:)   = 1.0;
gm(end,:) = 1.0;
% search across inner, *original* grid
for i = 2:imax+1
    for j = 2:jmax+1
        % test for cell being ocean
        if gm(j,i)
            cnt = gm(j+1,i) + gm(j-1,i) + gm(j,i+1) + gm(j,i-1);
            % test for cell surrounded by land (in 3 edge directions)
            if (cnt == 1)
                gc(j-1,i-1) = 1.0;
            end
        end
    end
end
grid_cells = gc;
%
% *********************************************************************** %
% *** END *************************************************************** %
% *********************************************************************** %
