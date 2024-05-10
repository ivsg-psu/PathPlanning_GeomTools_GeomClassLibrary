function agreement_indices = fcn_geometry_findAgreementOfPointsToCubicPoly(points, source_points, fittedParameters, transverse_tolerance, varargin)
%% fcn_geometry_findAgreementOfPointsToCubicPoly 
%
% Given a set of XY points, source points, and fitted parameters, finds the
% indicies of the points that are within a transverse_tolerance distance
% away from the cubic polynomial curve.
%
% FORMAT:
%
% agreement_indices = fcn_geometry_findAgreementOfPointsToCubicPoly(points, source_points, fittedParameters, transverse_tolerance, (fig_num))
%
% INPUTS:
%
%      points: a Nx2 vector where N is the number of points, but at least 2\
%      rows.
%
%      source_points: These are the points that are used to fit a cubic
%      polynomial curve using "polyfit"
%
%      fittedParameters: These are fitted parameters of the cubic
%      polynomial curve
%
%      transverse_tolerance: the orthogonal distance between the points and
%      the linear vector fit that indicate whether a point "belongs" to the
%      fit. A point belongs to the fit if the transverse distance is less
%      than or equal to the transverse_tolerance
%
%      (OPTIONAL INPUTS)
%
%      fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed.
%
% OUTPUTS:
%
%      agreement_indicies: the indicies of the points that are within
%      agreement of the best-fit parameters, given the transverse and
%      station tolerance settings. 
% 
% DEPENDENCIES:
%  
%      (NONE)
%
% EXAMPLES:
%
% See the script: script_test_fcn_geometry_findAgreementOfPointsToCubicPoly
% for a full test suite.

% This function was written on 2024_05_10 by Aneesh Batchu
% Questions or comments? abb6486@psu.edu or sbrennan@psu.edu

% Revision history:
% 2024_05_10 - Aneesh 
% -- wrote the code

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==5 && isequal(varargin{end},-1))
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
    debug_fig_num = 34838; %#ok<NASGU>
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

if 0==flag_max_speed
    if flag_check_inputs
        % Are there the right number of inputs?
        narginchk(5,6);

    end
end

% Does user want to specify fig_num?
fig_num = []; % Default is to have no figure
flag_do_plots = 0;
if (0==flag_max_speed) && (5<= nargin)
    temp = varargin{end};
    if ~isempty(temp)
        fig_num = temp;
        flag_do_plots = 1;
    end
end

%% Main Code starts from here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _       
%  |  \/  |     (_)      
%  | \  / | __ _ _ _ __  
%  | |\/| |/ _` | | '_ \ 
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find the slopes of cubic polynomial at each test source point
slopes_at_each_test_source_point = 3*fittedParameters(1,1)*source_points(:,1).^2 + 2*fittedParameters(1,2)*source_points(:,1) + fittedParameters(1,3);

% Find angle of inclination to find the unit_tangent_vector
theta = atan(slopes_at_each_test_source_point); 

% Unit tangent vectors of all test source points
unit_tangent_vectors = [sin(theta), cos(theta)]; 

% Unit orthogonal vectors of unit tangent vectors
unit_orthogonal_vectors = unit_tangent_vectors*[0 1; -1 0];

% Find the points on upper and lower boundaries
upper_boundary_points = source_points + unit_orthogonal_vectors.*transverse_tolerance;
lower_boundary_points = source_points - unit_orthogonal_vectors.*transverse_tolerance;

% Define the domain. These are the boundary vertices
polygon_vertices = [upper_boundary_points; lower_boundary_points(end:-1:1,:)];

% Determine which points lie inside the polygonal region (domain)
inlier_indices = inpolygon(points(:,1), points(:,2), polygon_vertices(:,1), polygon_vertices(:,2));
% outliers_indices = ~inliers_indices;

% Find the agreement indices based on inlier_indices
agreement_indices = find(inlier_indices==1);

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
    
    figure(fig_num)
    hold on;
    grid on;
    % axis equal
    xlabel('X [m]')
    ylabel('Y [m]')
    
    % plot the test points
    plot(points(:,1), points(:,2), 'r.', 'MarkerSize',30)

    % plot the source points of the cubic polynomial curve
    plot(source_points(:,1), source_points(:,2), 'b.', 'MarkerSize',30)

    % Plot the fitted polynomial
    x_fit = linspace(min(source_points(:,1)), max(source_points(:,1)), 100);
    y_fit = polyval(fittedParameters, x_fit);
    plot(x_fit, y_fit, 'r-', 'LineWidth', 2);
    
    % plot the upper boundary points
    plot(upper_boundary_points(:,1), upper_boundary_points(:,2), 'g.', 'MarkerSize',30)

    % plot the lower boundary points
    plot(lower_boundary_points(:,1), lower_boundary_points(:,2), 'm.', 'MarkerSize',30)

    % drawpolygon function is used to set the domain
    drawpolygon(gca,"Position",polygon_vertices);

    % Make axis slightly larger?
    if flag_rescale_axis
        temp = axis;
        %     temp = [min(points(:,1)) max(points(:,1)) min(points(:,2)) max(points(:,2))];
        axis_range_x = temp(2)-temp(1);
        axis_range_y = temp(4)-temp(3);
        percent_larger = 0.3;
        axis([temp(1)-percent_larger*axis_range_x, temp(2)+percent_larger*axis_range_x,  temp(3)-percent_larger*axis_range_y, temp(4)+percent_larger*axis_range_y]);
    end

end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file);
end

end

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


