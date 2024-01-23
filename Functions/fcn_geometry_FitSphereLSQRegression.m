%% fcn_LidarPoseEstimation_FitSphereLSQ
function [C_sphere,R_sphere,E_total] = fcn_geometry_FitSphereLSQRegression(XYZ_array)
% Fitting a sphere to points using least squares based on squared differences 
% of squared lengths and square radius
%
% FORMAT:
%
%       [C,R,E_total] = fcn_LidarPoseEstimation_FitSphereLSQ(XYZ_array)
%
% INPUTS:
%
%       XYZ_array: an Nx3 array contains X,Y and Z coordinates in Velodyne
%       Lidar Coordinate
%
% OUTPUTS:
%
%       C_sphere : a 1x3 vector contains the coordinate of the center of the sphere
%       R_sphere : the radius of the sphere
%       E_total  : ths sum of the squared differences of squared lengths 
%                  and square radius
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

flag_do_debug = 0; % Flag to plot the results for debugging
flag_check_inputs = 1; % Flag to perform input checking

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

if flag_check_inputs == 1
    % Are there the right number of inputs?
    narginchk(1,1);

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
dist_vec = vecnorm(diff_vec,2,2);
R_squre_sum = sum(dist_vec.^2);
R_squre = R_squre_sum/N_points;
R_sphere = sqrt(R_squre);
E_total = sum((dist_vec.^2 - R_sphere^2).^2);


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
% % Plot the inputs?    
% if flag_do_plots
% 
%     
%     % Nothing to do here!
% end

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
