%% fcn_LidarPoseEstimation_FitSphereLSQ
function [C_sphere,R_sphere,E_total, errors] = fcn_geometry_FitSphereLSQRegression(XYZ_array, varargin)
% Fitting a sphere to points using least squares based on squared differences 
% of squared lengths and square radius
%
% FORMAT:
%
%       [C,R,E_total,errors] = fcn_LidarPoseEstimation_FitSphereLSQ(XYZ_array)
%
% INPUTS:
%
%       XYZ_array: an Nx3 array contains X,Y and Z coordinates in Velodyne
%       Lidar Coordinate, with N>=3.
%
% OUTPUTS:
%
%       C_sphere : a 1x3 vector contains the coordinate of the center of the sphere
%       R_sphere : the radius of the sphere
%       E_total  : the sum of the squared differences of squared lengths 
%                  and square radius
%       errors   : the differences of point distances to center to fitted
%       radius
%
% DEPENDENCIES:
%
%      (none)
%
% EXAMPLES:
%
%       See the script:
%       script_test_fcn_LidarPoseEstimation_FitSphereLSQ.m for a full
%       test suite.
%
% This function was written on 2023_10_22 by X. Cao
% Questions or comments? xfc5113@psu.edu


% Revision history:
% 2023_10_22 by X. Cao
% -- start writing function
% 2024_01_23 by S. Brennan
% -- fixed header area (debugging)
% -- added varargin input for figures
% -- added input checking
% -- added plotting
% -- added errors output
% 2024_01_25 by S. Brennan
% -- added color-coded errors

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
    debug_fig_num = 34838;
else
    debug_fig_num = [];  
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

        % Check the XYZ_array input to have 3 columns, 3 or more rows of data
        fcn_DebugTools_checkInputsToFunctions(...
            XYZ_array, '3column_of_numbers',[3 4]);

    end
end

% Does user want to specify fig_num?
fig_num = []; % Default is to have no figure
flag_do_plots = 0;
if (0==flag_max_speed) && (2<= nargin)
    temp = varargin{end};
    if ~isempty(temp)
        fig_num = temp;
        flag_do_plots = 1;
    end
end

% Setup figures if there is debugging
if flag_do_debug
    fig_debug = 9999; 
else
    fig_debug = []; %#ok<*NASGU>
end

%% Write main code for plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _
%  |  \/  |     (_)
%  | \  / | __ _ _ _ __
%  | |\/| |/ _` | | '_ \
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N_points = size(XYZ_array,1);
XYZ_ave = mean(XYZ_array,1); % Compute the average of the data points

% Compute the covariance matrix P of the XYZ_diff = XYZ_array-XYZ_ave 
% and the right-hand side Q of the linear system P*(C_sphere - XYZ_ave) = Q.
XYZ_diff = XYZ_array - XYZ_ave; % Calculate the difference between each data point and the average of the data points
P_sum = zeros(3,3); % The sum of the covariance matrix P
Q_sum = zeros(3,1); % The sum of Q
P = XYZ_diff.'*XYZ_diff;
Q_1 = sum((XYZ_diff.').^2,1);
Q = sum(sum((XYZ_diff.').^2,1).*XYZ_diff.',2);

C_sphere = XYZ_ave + 1/2*(P\Q).';
diff_vec = C_sphere - XYZ_array;
dist_vec = sum(diff_vec.^2,2).^0.5;
R_squre_sum = sum(dist_vec.^2);
R_squre = R_squre_sum/N_points;
R_sphere = sqrt(R_squre);
errors = dist_vec - R_sphere;
E_total = sum(errors.^2);

%% Any debugging?
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
% Plot the inputs?    
if flag_do_plots
     % set up the figure and check if it has been used before. If not used
    % already, note that the axes will be re-scaled
    temp_h = figure(fig_num);
    flag_rescale_axis = 0;
    if isempty(get(temp_h,'Children'))
        flag_rescale_axis = 1;
        hold on;
        axis equal
        grid minor;
        view(3);

        xlabel('X [m]');
        ylabel('Y [m]');
        zlabel('Z [m]');

    else
        hold on;
        axis equal
    end

    
    % Plot the input points
    flag_show_error_in_color = 1;
    if flag_show_error_in_color==0
        plot3(XYZ_array(:,1), XYZ_array(:,2),  XYZ_array(:,3), 'k.','MarkerSize',20);
    else
        % This is from temp = colormap('turbo')
        error_colorMap = [   
            0.1900    0.0718    0.2322
            0.2737    0.3835    0.8449
            0.2093    0.6654    0.9760
            0.1034    0.8960    0.7150
            0.4483    0.9959    0.3694
            0.7824    0.9376    0.2033
            0.9800    0.7300    0.2216
            0.9643    0.4186    0.0964
            0.7941    0.1660    0.0143
            0.4796    0.0158    0.0106];
        max_error = max(abs(errors));
        error_indicies = max(1,round(abs(errors)./max_error * length(error_colorMap(:,1))));
        
        for ith_color = 1:length(error_colorMap(:,1))
            current_color = error_colorMap(ith_color,:);
            points_to_plot = XYZ_array(error_indicies==ith_color,:);
            plot3(points_to_plot(:,1), points_to_plot(:,2),  points_to_plot(:,3), '.','MarkerSize',20,'Color',current_color);
        end

    end

    % Plot the fitted sphere
    plot3(C_sphere(1,1), C_sphere(1,2), C_sphere(1,3), 'r+','MarkerSize',30);

    % Plot the sphere
    fcn_geometry_plotSphere(C_sphere, R_sphere, [0 0 1], fig_num);

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

end % Ends main function for fcn_LidarPoseEstimation_FitSphereLSQ

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
