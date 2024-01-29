function [test_points, true_base_points, true_projection_vectors, true_distances] = fcn_geometry_fillLineTestPoints(seed_points, M, sigma, varargin)
% fcn_geometry_fillLineTestPoints
% given N points, with N>=2, creates a set of M points per unit distance
% between these points randomly distributed with variance sigma.
%
% [test_points] = fcn_geometry_fillLineTestPoints(seed_points, M, sigma)
%
% INPUTS:
%
%      seed_points: a Nx2 or Nx3 vector where N is the number of points,
%      but at least 2.
%
%      M: the number of test points to generate per unit
%      distance.
%
%      sigma: athe standard deviation in points
%
% OUTPUTS:
%
%      test_points: a list of test points used to test regression fitting
%
%      true_base_points, true_projection_vectors, true_distances: the true
%      values of the fitting
%
% DEPENDENCIES:
%
%      fcn_DebugTools_checkInputsToFunctions
%      fcn_geometry_calcUnitVector
%      fcn_geometry_calcOrthogonalVector
%
% EXAMPLES:
%      
% See the script: script_test_fcn_geometry_fillLineTestPoints
% for a full test suite.
%
% NOTE: This function does NOT work for fitting all lines.
%
% This function was written on 2023_12_05 by S. Brennan
% Questions or comments? sbrennan@psu.edu 

% Revision history:
% 2023_12_05 
% -- wrote the code
% 2024_01_03 - S. Brennan
% -- added fast mode option
% -- added environmental variable options
% 2024_01_23 - S. Brennan
% -- added support for 3 dimensional points

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==4 && isequal(varargin{end},-1))
    flag_do_debug = 0; % Flag to plot the results for debugging
    flag_check_inputs = 0; % Flag to perform input checking
    flag_max_speed = 1;
else
    % Check to see if we are externally setting debug mode to be "on"
    flag_do_debug = 0; % Flag to plot the results for debugging
    flag_check_inputs = 1; % Flag to perform input checking
    MATLABFLAG_GEOMETRY_FLAG_CHECK_INPUTS = getenv("MATLABFLAG_GEOMETRY_FLAG_CHECK_INPUTS");
    MATLABFLAG_GEOMETRY_FLAG_DO_DEBUG = getenv("MATLABFLAG_GEOMETRY_FLAG_DO_DEBUG");
    if ~isempty(MATLABFLAG_GEOMETRY_FLAG_CHECK_INPUTS) && ~isempty(MATLABFLAG_GEOMETRY_FLAG_DO_DEBUG)
        flag_do_debug = str2double(MATLABFLAG_GEOMETRY_FLAG_DO_DEBUG);
        flag_check_inputs  = str2double(MATLABFLAG_GEOMETRY_FLAG_CHECK_INPUTS);
    end
end

if flag_do_debug
    st = dbstack; %#ok<*UNRCH>
    fprintf(1,'STARTING function: %s, in file: %s\n',st(1).name,st(1).file);
    debug_fig_num = 34838;
else
    debug_fig_num = []; %#ok<NASGU>
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

if (0==flag_max_speed)
    if flag_check_inputs == 1
        % Are there the right number of inputs?
        narginchk(3,4);

        % Check the points input to be length greater than or equal to 2
        fcn_DebugTools_checkInputsToFunctions(...
            seed_points, '2or3column_of_numbers',[2 3]);

    end
end

% Does user want to show the plots?
flag_do_plots = 0;
if (0==flag_max_speed)
    if (4 == nargin)
        temp = varargin{end};
        if ~isempty(temp)
            fig_num = temp;
            figure(fig_num);
            flag_do_plots = 1;
        end
    elseif flag_do_debug
        fig = figure;
        fig_num = fig.Number;
    end
end

%% Solve for the circle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _       
%  |  \/  |     (_)      
%  | \  / | __ _ _ _ __  
%  | |\/| |/ _` | | '_ \ 
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_segments = length(seed_points(:,1)) -1;

distances = sum((seed_points(2:end,:) - seed_points(1:end-1,:)).^2,2).^0.5;
unit_vectors = fcn_geometry_calcUnitVector(seed_points(2:end,:)-seed_points(1:end-1,:),-1);
unit_orthogonals = fcn_geometry_calcOrthogonalVector(unit_vectors, [], -1);

test_points = [];
for ith_point = 1:N_segments
    projection_distances = (0:(1/M):distances(ith_point))';
    N_points = length(projection_distances);
    orthogonal_distances = randn(N_points,1)*sigma;

    if length(seed_points(1,:))==2        
        test_points = [test_points; seed_points(ith_point,:) + projection_distances*unit_vectors(ith_point,:) + orthogonal_distances.*unit_orthogonals(ith_point,:)]; %#ok<AGROW>
    else
        segment_orthogonals = fcn_geometry_calcOrthogonalVector(ones(N_points,1)*unit_vectors(ith_point,:),[], -1);
        test_points = [test_points; seed_points(ith_point,:) + projection_distances*unit_vectors(ith_point,:) + orthogonal_distances.*segment_orthogonals]; %#ok<AGROW>
    end
end

true_base_points = seed_points(1:end-1,:);
true_projection_vectors = unit_vectors;
true_distances = distances;

%% Plot the results (for debugging)?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____       _                 
%  |  __ \     | |                
%  | |  | | ___| |__  _   _  __ _ 
%  | |  | |/ _ \ '_ \| | | |/ _` |
%  | |__| |  __/ |_) | |_| | (_| |
%  |_____/ \___|_.__/ \__,_|\__, |
%                            __/ |
%                           |___/ 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flag_do_plots


    temp_h = figure(fig_num);
    flag_rescale_axis = 0;
    if isempty(get(temp_h,'Children'))
        flag_rescale_axis = 1;
    end        

    hold on;
    grid on;
    axis equal
    xlabel('X [m]')
    ylabel('Y [m]')

    % Plot the input vectors alongside the unit vectors
    if length(seed_points(1,:))==2 % 2D vectors
        % Plot the input points
        plot(seed_points(:,1),seed_points(:,2),'r.-','MarkerSize',20);

        % Plot the results
        plot(test_points(:,1),test_points(:,2),'b.','MarkerSize',10);


    else % 3D vectors
        zlabel('Z [m]')

        view(3);
        % Plot the input points
        plot3(seed_points(:,1),seed_points(:,2),seed_points(:,3),'r.-','MarkerSize',20);

        % Plot the results
        plot3(test_points(:,1),test_points(:,2),test_points(:,3),'b.','MarkerSize',10);

    end

    % Make axis slightly larger?
    if flag_rescale_axis
        temp = axis;
        %     temp = [min(points(:,1)) max(points(:,1)) min(points(:,2)) max(points(:,2))];
        axis_range_x = temp(2)-temp(1);
        axis_range_y = temp(4)-temp(3);
        percent_larger = 0.3;
        axis([temp(1)-percent_larger*axis_range_x, temp(2)+percent_larger*axis_range_x,  temp(3)-percent_larger*axis_range_y, temp(4)+percent_larger*axis_range_y]);
    end

    
end % Ends check if plotting

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file);
end

end % Ends main function

%% Functions follow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   ______                _   _
%  |  ____|              | | (_)
%  | |__ _   _ _ __   ___| |_ _  ___  _ __  ___
%  |  __| | | | '_ \ / __| __| |/ _ \| '_ \/ __|
%  | |  | |_| | | | | (__| |_| | (_) | | | \__ \
%  |_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|___/
%
% See: https://patorjk.com/software/taag/#p=display&f=Big&t=Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%§







