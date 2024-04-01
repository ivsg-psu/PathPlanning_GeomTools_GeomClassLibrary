function [radius, arcCenter, arcLength, radial_fitting_error] = fcn_geometry_fitArcConstrainedStart(input_points, varargin)
%% fcn_geometry_fitArcConstrainedStart
% Given a portion of an arc where the arc must start at a particular
% orientation, performs constrained least-error regression minimizing
% radial error to find the best-fit radius and center of the arc. Also
% calculates the radial fitting error for each point.
% 
% Format: 
% [radius, arcCenter, arcLength, radial_fitting_error] = fcn_geometry_fitArcConstrainedStart(input_points, (initial_rotation), (initial_offset), (fig_num))
%
% INPUTS:
%      input_points: an [Nx2] matrix of N different [x y] points assumed to
%      be in sequence. Note: the function may break, particularly in the
%      calculation of the arcLength, if the points are not in sequence
%
%      (OPTIONAL INPUTS)
% 
%      initial_rotation: the initial angle of the arc's start point, in
%      radians. If nothing is entered, the arc is assumed to start with a
%      vertical orientation, e.g. initially going straight "up" the y-axis.
%
%      initial_offset: the initial position as [x y] of the arc's start
%      point. If nothing is entered, the arc is assumed to start at the
%      origin.
%
%      fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed.
%
% OUTPUTS:
%
%      radius: the radius of the arc fit
%
%      arcCenter: the [x y] location of the arc center
%
%      arcLength: the length of the arc, in radians
%
%      radial_fitting_error: an [Nx1] column of the radial fitting error,
%      where N is the number of input points.
%
% DEPENDENCIES:
%
%      (none)
%
% EXAMPLES:
%      
% See the script: script_test_fcn_geometry_fitArcConstrainedStart
% for a full test suite.
%
% This function was written on 2024_03_30 by S. Brennan
% Questions or comments? sbrennan@psu.edu 

% Revision history:
% 2024_03_30 - S Brennan
% -- wrote the code

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
    if flag_check_inputs == 1
        % Are there the right number of inputs?
        narginchk(1,4);

        % % Check the source_points input to be length exactly equal to 3
        % fcn_DebugTools_checkInputsToFunctions(...
        %     source_points, '2column_of_numbers',[3 3]);
        %
        % % Check the associated_points_in_domain input to be length greater
        % % than or equal to 3
        % fcn_DebugTools_checkInputsToFunctions(...
        %     associated_points_in_domain, '2column_of_numbers',[3 4]);
    end
end

% Does user want to specify initial_rotation?
initial_rotation = 0;
if (2<=nargin)
    temp = varargin{1};
    if ~isempty(temp)
        initial_rotation = temp;
    end
end

% Does user want to specify initial_rotation?
initial_offset = [0 0];
if (3<=nargin)
    temp = varargin{2};
    if ~isempty(temp)
        initial_offset = temp;
    end
end

% Does user want to specify fig_num?
fig_num = []; % Default is to have no figure
flag_do_plots = 0;
if  (0==flag_max_speed) && (4<= nargin)
    temp = varargin{end};
    if ~isempty(temp)
        fig_num = temp;
        flag_do_plots = 1;
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

input_points_translated = input_points - initial_offset;
input_points_translated_rotated = input_points_translated*[cos(initial_rotation) -sin(initial_rotation); sin(initial_rotation) cos(initial_rotation)];

% Perform regression fit - this was solved using a LOT of calculus (6 pages
% of derivations)
sum_y4 = sum(input_points_translated_rotated(:,2).^4,1);
sum_x2y2 = sum(input_points_translated_rotated(:,1).^2.*input_points_translated_rotated(:,2).^2,1);
sum_x1y2 = sum(input_points_translated_rotated(:,1).*input_points_translated_rotated(:,2).^2,1);
x0 = (sum_y4 + sum_x2y2)/(2*sum_x1y2);

% Save outputs
local_arcCenter = [x0 0];
radius = abs(x0);
radial_fitting_error = sum((input_points_translated_rotated-local_arcCenter).^2,2).^0.5 - radius;

% Transform back
arcCenter = local_arcCenter*[cos(initial_rotation) sin(initial_rotation); -sin(initial_rotation) cos(initial_rotation)] + initial_offset;

% Use test points to find the arc length
test_points = [input_points(1:2,:); input_points(end,:)];
angles = atan2((test_points(:,2) - arcCenter(1,2)),(test_points(:,1) - arcCenter(1,1)));
[~, arcLength, ~, ~, ~]  = fcn_geometry_arcAngleFrom3Points([cos(angles(1,1)) sin(angles(1,1))], [cos(angles(2,1)) sin(angles(2,1))], [cos(angles(end,1)) sin(angles(end,1))] );
arcLength = abs(arcLength);



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

    % Plot the results in point space
    temp_h = figure(fig_num);
    flag_rescale_axis = 0;
    if isempty(get(temp_h,'Children'))
        flag_rescale_axis = 1;
    end    

    % Get the color ordering?
    try
        color_ordering = orderedcolors('gem12');
    catch
        color_ordering = colororder;
    end

    N_colors = length(color_ordering(:,1));

    hold on;
    grid on;
    axis equal;

    % Plot the fits    
    ith_domain = 1;
    current_color = color_ordering(mod(ith_domain,N_colors)+1,:); 
      
    
    % Plot the associated_points_in_domain
    plot(input_points(:,1),input_points(:,2),'.','MarkerSize',20,'Color',current_color);

    % Plot the fit
    plot(arcCenter(:,1),arcCenter(:,2),'+','Color','b','MarkerSize',40);

    angles = atan2((input_points(:,2) - arcCenter(1,2)),(input_points(:,1) - arcCenter(1,1)));
    [~, plotting_arcLength, ~, ~, start_angles_in_radians]  = fcn_geometry_arcAngleFrom3Points([cos(angles(1,1)) sin(angles(1,1))], [cos(angles(2,1)) sin(angles(2,1))], [cos(angles(end,1)) sin(angles(end,1))] );

    max_angle = start_angles_in_radians+plotting_arcLength;

    if plotting_arcLength>0
        angle_sweep = [(start_angles_in_radians:1*pi/180:(start_angles_in_radians+plotting_arcLength))'; max_angle];
    else
        angle_sweep = [(start_angles_in_radians:-1*pi/180:(start_angles_in_radians+plotting_arcLength))'; max_angle];
    end
    fitted_points = abs(x0)*[cos(angle_sweep) sin(angle_sweep)] + arcCenter;
    plot(fitted_points(:,1),fitted_points(:,2),'-','MarkerSize',40,'Color',[0 0 1],'LineWidth',3);

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


