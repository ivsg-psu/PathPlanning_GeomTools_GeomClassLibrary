function [gridIndices,gridDomains,gridCenters] = fcn_geometry_separatePointsIntoGrids(inputPoints, gridSize, gridBoundaries)
% Introduction:
% This function aims to seprate given points in X/XY/XYZ format into a
% user defined grid size and grid boundaries and output a cell array
% containing indices of points in the domain,the domain in AABB format and
% the grid centers
% Example: a set of XYZ points, grid space of two meters, the function will
% seprate all the points in [0,2m], [2m,4m], etc. If the points are in 2D
% there should be an option in the function to create the grid in x-axis or
% in y-axis
%
% FORMAT:
%
%       [gridIndices,gridDomains,gridCenters] = fcn_geometry_separatePointsIntoGrids(inputPoints, gridSize, gridBoundaries)
%
% INPUTS:
%
%      (OPTIONAL INPUTS)
%
%     none
%
% OUTPUTS:
%
%     gridIndices : a cell array containing indices of points for specific
%     domains
%
%     gridDomains : matrix in the format of AABB (i.e axis aligned bounding
%     box format) that gives the domains corresponding to the indices of
%     points present in the domain
%
%     gridCenters : matrix in 2D or 3D that gives the respective
%     gridDomain centers
%
%
% DEPENDENCIES:
%
%      (none)
%
% EXAMPLES:
%
%       See the script:
%       scripts_test_fcn_geometry_separatePointsIntoGrids.m for a full
%       test suite.
%
% This function was written on 2024_01_22 by V. Wagh
% Questions or comments? vbw5054@psu.ed


% Revision history:
% 2024_01_22 by V. Wagh (vbw5054@psu.edu)
% -- start writing code
% 2024_01_24 by V. Wagh
% -- changed inputs to include gridBoundaries
% -- changed outputs

flag_do_debug = 0; % Flag to plot the results for debugging
flag_check_inputs = 1; % Flag to perform input checking

if flag_do_debug
    st = dbstack; %#ok<*UNRCH>
    fprintf(1,'STARTING function: %s, in file: %s\n',st(1).name,st(1).file);
    debug_fig_num = 34838;
else
    debug_fig_num = [];  % no plotting yet
end
%% check input arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____                   _
%  |_   _|                 | |
%    | |  _ __  _ __  _   _| |_ ___
%    | | | '_ \| '_ \| | | | __/ __|
%   _| |_| | | | |_) | |_| | |_\__ \
%  |_____|_| |_| .__/ \__,_|\__|___/
%              | |
%              |_|
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if flag_check_inputs == 1
    % Are there the right number of inputs?
    narginchk(3,3); % no optional inputs
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _
%  |  \/  |     (_)
%  | \  / | __ _ _ _ __
%  | |\/| |/ _` | | '_ \
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Steps
% 1. Take inputs: array of random points to be seprated in 1D 2D or 3D
% format
% 2. Grid size in  meters
% 3. Are points in 2D or 3D
% 4. If pints are in 2D, do you want grid in x-axis or y-axis?

% Calculate the grid boundaries from the gridBoundaries input
xMin = gridBoundaries(1,1);
xMax = gridBoundaries(1,2);
if size(gridBoundaries,1) > 1
    yMin = gridBoundaries(2,1);
    yMax = gridBoundaries(2,2);

    if size(gridBoundaries,1) == 3
        zMin = gridBoundaries(3,1);
        zMax = gridBoundaries(3,2);
    else
        zMin = 0;
        zMax = 0;
    end
else
    yMin = 0;
    yMax = 0;
    zMin = 0;
    zMax = 0;
end

% % what is the lowest number betwween xMax, yMax, zMax
% axesMin = max(xMin,yMin,zMin);
%
% % what is the highest number betwween xMax, yMax, zMax
% axesMax = max(xMax,yMax,zMax);


% Calculate the number of grids in each dimension
numGridsX = ceil((xMax - xMin) / gridSize);
numGridsY = ceil((yMax - yMin) / gridSize);
numGridsZ = ceil((zMax - zMin) / gridSize);

% total number of domains
if size(gridBoundaries,1) == 1
    totalDomains = numGridsX;
elseif size(gridBoundaries,1) == 2
    totalDomains = numGridsX * numGridsY;
elseif size(gridBoundaries,1) == 3
    totalDomains = numGridsX * numGridsY * numGridsZ;
end


% % once we have the gridBoundaries in x, y and z we can go ahead and
% % calculate the domains in the AABB format
%
% % iterating through the grids

count_var = 0;
if size(gridBoundaries,1) == 1
    for i_1D = 1:numGridsX
        % for 1D (only x is evaluated)
        % Define the boundaries of the current grid
        gridXMin = xMin + (i_1D - 1) * gridSize;
        gridXMax = xMin + i_1D * gridSize;
        gridDomains{count_var + i_1D} = [gridXMin gridXMax];

        % Find indices of points within the current grid
        pointsInGrid = find(inputPoints(:, 1) >= gridXMin & inputPoints(:, 1) < gridXMax);

        gridIndices{count_var + i_1D} = pointsInGrid;
        gridCenters{count_var + i_1D} = [mean([gridXMin, gridXMax])];
    end

    % for 2D (only XY is evaluated)
elseif size(gridBoundaries,1) == 2
    for i_1D = 1:numGridsX
        for j_2D = 1:numGridsY
            gridXMin = xMin + (i_1D - 1) * gridSize;
            gridXMax = xMin + i_1D * gridSize;
            gridYMin = yMin + (j_2D - 1) * gridSize;
            gridYMax = yMin + j_2D * gridSize;
            gridDomains{count_var + j_2D} = [gridXMin gridXMax gridYMin gridYMax];

            % Find indices of points within the current grid
            pointsInGrid = find(inputPoints(:, 1) >= gridXMin & inputPoints(:, 1) < gridXMax & ...
                inputPoints(:, 2) >= gridYMin & inputPoints(:, 2) < gridYMax);

            gridIndices{count_var + j_2D} = pointsInGrid;
            gridCenters{count_var + j_2D} = [mean([gridXMin, gridXMax]), mean([gridYMin, gridYMax])];
        end
        count_var = count_var + j_2D;
    end

    % for 3D (XYZ is evaluated)
elseif size(gridBoundaries,1) == 3
    for i_1D = 1:numGridsX
        for j_2D = 1:numGridsY
            for k_3D = 1:numGridsZ
                gridXMin = xMin + (i_1D - 1) * gridSize;
                gridXMax = xMin + i_1D * gridSize;
                gridYMin = yMin + (j_2D - 1) * gridSize;
                gridYMax = yMin + j_2D * gridSize;
                gridZMin = zMin + (k_3D - 1) * gridSize;
                gridZMax = zMin + k_3D * gridSize;
                gridDomains{count_var + j_2D + k_3D} = [gridXMin gridXMax gridYMin gridYMax gridZMin gridZMax];

                % Find indices of points within the current grid
                pointsInGrid = find(inputPoints(:, 1) >= gridXMin & inputPoints(:, 1) < gridXMax & ...
                    inputPoints(:, 2) >= gridYMin & inputPoints(:, 2) < gridYMax & ...
                    inputPoints(:, 2) >= gridZMin & inputPoints(:, 2) < gridZMax);

                gridIndices{count_var + j_2D + k_3D} = pointsInGrid;
                gridCenters{count_var + j_2D + k_3D} = [mean([gridXMin, gridXMax]), mean([gridYMin, gridYMax]), mean([gridZMin, gridZMax])];
            end
            count_var = count_var + j_2D;
        end
        count_var = count_var + j_2D + k_3D;
    end
end
end


%% old code
%     % Initialize cell array to store grid data
%     gridData = cell(numGridsX, numGridsY);
%
%     % Iterate over each grid
%     for i = 1:numGridsX
%         for j = 1:numGridsY
%             % Define the boundaries of the current grid
%             gridXMin = xMin + (i - 1) * gridSize;
%             gridXMax = xMin + i * gridSize;
%             gridYMin = yMin + (j - 1) * gridSize;
%             gridYMax = yMin + j * gridSize;
%
%             % Find indices of points within the current grid
%             pointsInGrid = find(inputPoints(:, 1) >= gridXMin & inputPoints(:, 1) < gridXMax & ...
%                                 inputPoints(:, 2) >= gridYMin & inputPoints(:, 2) < gridYMax);
%
%             % Store grid data in the cell array
%             gridData{i, j}.indices = pointsInGrid;
%             gridData{i, j}.boundaries = [gridXMin, gridXMax, gridYMin, gridYMax];
%             gridData{i, j}.center = [mean([gridXMin, gridXMax]), mean([gridYMin, gridYMax])];
%         end
%     end
% end
