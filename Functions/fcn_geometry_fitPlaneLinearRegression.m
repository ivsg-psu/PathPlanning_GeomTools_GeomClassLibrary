function [parameters, standard_deviation_in_z, z_fit] = ...
    fcn_geometry_fitPlaneLinearRegression(points,varargin)
% fcn_geometry_fitPlaneLinearRegression
% Finds the coefficients, A, B, and C in the equation:
%
% z = Ax + By + C
%
% FORMAT: 
%
% [parameters, standard_deviation_in_z, z_fit] = ...
%    fcn_geometry_fitPlaneLinearRegression(points,(fig_num))
%
% INPUTS:
%
%      points: a Nx3 vector where N is the number of points, length N>=3. 
%
%      (OPTIONAL INPUTS)
% 
%      fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed.
%
% OUTPUTS:
%
%      params: the [1 x 3] matrix of the parameters [A B C]
%
%      standard_deviation_in_z: the standard deviation in the z-error of
%      the 
%  
%     z_fit: the model-fit z values
%
% DEPENDENCIES:
%
%      fcn_DebugTools_checkInputsToFunctions
%
% EXAMPLES:
%       
% See the script: script_test_fcn_geometry_fitPlaneLinearRegression
% for a full test suite.
%
% This function was written on 2024_01_19 by S. Brennan
% Questions or comments? sbrennan@psu.edu 


% Revision history:
% 2024_01_19 - S. Brennan
% -- wrote the code

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
 
% The code below uses the least-y-squared solution of:
% z = Ax + By + C
%
% z = [x y 1]*[A B C]'
%
% Let X = [x y 1], P = [A B C]'
%
% then
%   z = X*P
% 
% which, solved for P:
%
% inv(X'*X)*(X'z) = P

x = points(:,1);
y = points(:,2);
N_points = length(x(:,1));

X = [x, y, ones(length(points(:,1)),1)];
z = points(:,3);

P = (X'*X)\(X'*z);

parameters = P;

z_fit = [x y ones(N_points,1)]*parameters; % Solve for z vertices data
z_error = z - z_fit;
standard_deviation_in_z = std(z_error, 0, 1);


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


    % Plot the results in point space
    temp_h = figure(fig_num);
    flag_rescale_axis = 0;
    if isempty(get(temp_h,'Children'))
        flag_rescale_axis = 1;
    end      

    hold on;
    grid on;
    axis equal;

    % Plot the points
    plot3(points(:,1),points(:,2),points(:,3),'k.','MarkerSize',20);
    plot3(points(:,1),points(:,2),z_fit,'.','MarkerSize',20);
    view(3);

    % Make axis slightly larger?
    if flag_rescale_axis
        temp_axis = axis;
        %     temp = [min(points(:,1)) max(points(:,1)) min(points(:,2)) max(points(:,2))];
        axis_range_x = temp_axis(2)-temp_axis(1);
        axis_range_y = temp_axis(4)-temp_axis(3);
        percent_larger = 0.3;
        axis([temp_axis(1)-percent_larger*axis_range_x, temp_axis(2)+percent_larger*axis_range_x,  temp_axis(3)-percent_larger*axis_range_y, temp_axis(4)+percent_larger*axis_range_y]);
    else
        temp_axis = axis;
    end


    % Plot the plane
    % x = [1 -1 -1 1]; % Generate data for x vertices
    x = [temp_axis(2) temp_axis(1) temp_axis(1) temp_axis(2)]';

    % y = [1 1 -1 -1]; % Generate data for y vertices
    y = [temp_axis(4) temp_axis(4) temp_axis(3) temp_axis(3)]';

    N_points = length(x(:,1));
    z = [x y ones(N_points,1)]*parameters; % Solve for z vertices data

    h_patch = patch(x, y, z, [0 0 1],'FaceAlpha',0.1); %#ok<NASGU>


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

