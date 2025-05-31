function [vector_root, unit_vector, standard_deviation_in_plane_orthogonals, plane_distances] = ...
    fcn_geometry_fitVectorToPlane(points,varargin)
% fcn_geometry_fitVectorToPlane
% Finds the vector that is normal to a set of points that define a 3D
% plane. The method is to find the vector fitting XY, YZ, and XZ data, and
% then using cross products between the two largest vectors to find the
% normal.
%
% FORMAT: 
%
% [vector_root, unit_vector] = fcn_geometry_fitVectorToPlane(points)
%
% INPUTS:
%
%      points: a Nx3 vector where N is the number of points, length N>=2. 
%
%      (OPTIONAL INPUTS)
% 
%      fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed.
%
% OUTPUTS:
%
%      vector_root: the [1 x 3] matrix of the (x,y) point representing the
%      root of the vector, wherein the root is the point on the vector
%      closest to the origin.
%
%      unit_vector: the [1 x 3] matrix of the (deltax,deltay,deltaz) length of the
%      unit vector attached to the root.
%
%     standard_deviation_in_plane_orthogonals: the standard deviation in
%     the point fitting error in the direction of the unit_vector
%
%     plane_distances: the distances of each of the N points to the plane,
%     measured orthogonally from the plane. Returned as an [N x 1 ] vector.
%
% DEPENDENCIES:
%
%      fcn_DebugTools_checkInputsToFunctions
%      fcn_geometry_fitVectorToNPoints
%
% EXAMPLES:
%      
%      % BASIC example
%      points = [2 3 4; 4 5 6; 1 1 1];
%      [vector_root, unit_vector] = fcn_geometry_fitVectorToPlane(points)
% 
% See the script: script_test_fcn_geometry_fitVectorToPlane
% for a full test suite.
%
% This function was written on 2024_01_27 by S. Brennan
% Questions or comments? sbrennan@psu.edu 


% Revision history:
% 2024_01_27 
% -- wrote the code
% -- modified from fcn_geometry_fitVectorToNPoints


%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==2 && isequal(varargin{end},-1))
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
    if flag_check_inputs == 1
        % Are there the right number of inputs?
        narginchk(1,2);

        % Check the points input to be length greater than or equal to 2
        fcn_DebugTools_checkInputsToFunctions(...
            points, '3column_of_numbers',[3 4]);

    end
end

% Does user want to show the plots?
flag_do_plot = 0;
if (0==flag_max_speed) && (2 == nargin) 
    temp_axis = varargin{1};
    if ~isempty(temp_axis)
        fig_num = temp_axis;
        figure(fig_num);
        flag_do_plot = 1;
    end
else
    if flag_do_debug
        fig = figure; 
        fig_num = fig.Number;
        flag_do_plot = 1;
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
x = points(:,1);
y = points(:,2);
z = points(:,3);

fitting_fig_num = -1; % Force speed mode

% Find the XYZ vector solution
[~, standard_deviation_in_z_XYZ, ~, unit_vector_XYZ, base_point_XYZ, standard_deviation_in_plane_orthogonals_XYZ, plane_distances_XYZ] = fcn_geometry_fitPlaneLinearRegression([x y z],fitting_fig_num);

% Find the XZY vector
[~, standard_deviation_in_z_XZY, ~, unit_vector_XZY, base_point_XZY, standard_deviation_in_plane_orthogonals_XZY, plane_distances_XZY] = fcn_geometry_fitPlaneLinearRegression([x z y],fitting_fig_num);

% Find the YZX vector
[~, standard_deviation_in_z_YZX, ~, unit_vector_YZX, base_point_YZX, standard_deviation_in_plane_orthogonals_YZX, plane_distances_YZX] = fcn_geometry_fitPlaneLinearRegression([y z x],fitting_fig_num);

% Perform reassembly
unit_vector_XYZ_corrected = unit_vector_XYZ;
unit_vector_XZY_corrected = [unit_vector_XZY(1) unit_vector_XZY(3) unit_vector_XZY(2)] ;
unit_vector_YZX_corrected = [unit_vector_YZX(3) unit_vector_YZX(1) unit_vector_YZX(2)];

base_point_XYZ_corrected = base_point_XYZ;
warning('base point calculated in fcn_geometry_fitVectorToPlane may be in error, as plane fit regression function has been updated.')
base_point_XZY_corrected = [base_point_XZY(1) base_point_XZY(3) base_point_XZY(2)] ;
base_point_YZX_corrected = [base_point_YZX(3) base_point_YZX(1) base_point_YZX(2)];

[~,min_std_z_index] = min([standard_deviation_in_z_XYZ, standard_deviation_in_z_XZY, standard_deviation_in_z_YZX]); %#ok<ASGLU>
[~,min_std_ortho_index] = min([standard_deviation_in_plane_orthogonals_XYZ, standard_deviation_in_plane_orthogonals_XZY, standard_deviation_in_plane_orthogonals_YZX]);

% if min_std_z_index~=min_std_ortho_index
%     warning('on','backtrace');
%     warning('minimum direction for standard deviations in z does not match minimum directions for standard deviations in plane orthogonals.');
% end

% Find the root and unit vector depending on which situation gives minimum
% standard deviation
switch min_std_ortho_index
    case 1
        vector_root = base_point_XYZ_corrected;
        unit_vector = unit_vector_XYZ_corrected;
        standard_deviation_in_plane_orthogonals = standard_deviation_in_plane_orthogonals_XYZ;
        plane_distances = plane_distances_XYZ;
    case 2
        vector_root = base_point_XZY_corrected;
        unit_vector = unit_vector_XZY_corrected;
        standard_deviation_in_plane_orthogonals = standard_deviation_in_plane_orthogonals_XZY;
        plane_distances = plane_distances_XZY;
    case 3
        vector_root = base_point_YZX_corrected;
        unit_vector = unit_vector_YZX_corrected;
        standard_deviation_in_plane_orthogonals = standard_deviation_in_plane_orthogonals_YZX;
        plane_distances = plane_distances_YZX;
    otherwise 
        warning('on','backtrace');
        warning('minimum direction for standard deviations in z does not match minimum directions for standard deviations in plane orthogonals.');
        error('unexpected situation occurred, unable to find minimum standard deviations.');
end


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
if flag_do_plot
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
    ylabel('Z [m]')
    view(3)

    % Plot the input points    
    plot3(points(:,1),points(:,2),points(:,3),'k.');
    
    % Plot the root
    plot3(vector_root(:,1),vector_root(:,2),vector_root(:,3),'b.','Markersize',50);
    
    % Plot the unit vector
    quiver3(vector_root(1,1),vector_root(1,2),vector_root(1,3), unit_vector(1,1),unit_vector(1,2),unit_vector(1,3),0,'b','Linewidth',3);

    % Make axis slightly larger?
    temp_axis = axis;
    if flag_rescale_axis
        %     temp = [min(points(:,1)) max(points(:,1)) min(points(:,2)) max(points(:,2))];
        axis_range_x = temp_axis(2)-temp_axis(1);
        axis_range_y = temp_axis(4)-temp_axis(3);
        percent_larger = 0.3;
        axis([temp_axis(1)-percent_larger*axis_range_x, temp_axis(2)+percent_larger*axis_range_x,  temp_axis(3)-percent_larger*axis_range_y, temp_axis(4)+percent_larger*axis_range_y]);
    end


    % Plot the plane
    up_direction = fcn_geometry_calcOrthogonalVector(unit_vector);
    right_direction = cross(unit_vector,up_direction);

    distances_up    = sum(up_direction.*points,2);
    distances_right = sum(right_direction.*points,2);

    max_up = max(distances_up);
    min_up = min(distances_up);
    max_right = max(distances_right);
    min_right = min(distances_right);

    plane_boundaries = ones(4,1)*vector_root + [
        min_up*up_direction+min_right*right_direction;
        min_up*up_direction+max_right*right_direction;
        max_up*up_direction+max_right*right_direction;
        max_up*up_direction+min_right*right_direction;
        ];
    h_patch = patch(plane_boundaries(:,1), plane_boundaries(:,2), plane_boundaries(:,3), [0 0 1],'FaceAlpha',0.1); %#ok<NASGU>


    
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


%% OLD METHOD (WRONG)

% % Find the mean values of the differences in the points, row-wise
% diff_points = diff(points);
% if length(diff_points(:,1))>1
%     mean_diff_points = mean(diff_points);
% else
%     mean_diff_points = diff_points;
% end
% 
% % Calculate x1, x2, y1, y2 by extracting each from points matrix
% x1 = points(1:end-1,1);
% y1 = points(1:end-1,2);
% x2 = points(2:end,1);
% y2 = points(2:end,2);
% 
% % Grab A, B, C values
% A = -1*mean_diff_points(1,2);
% B =    mean_diff_points(1,1);
% C = mean(x1.*y2 - y1.*x2);

