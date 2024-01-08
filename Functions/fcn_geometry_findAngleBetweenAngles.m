function [isAngleBetween]  = fcn_geometry_findAngleBetweenAngles(start_angle_in_radians, end_angle_in_radians, direction, angles_to_test_in_radians, varargin)
% fcn_geometry_findAngleBetweenAngles checks whether angles lie between two
% different angles
%
% FORMAT:
%
% [isAngleBetween]  = fcn_geometry_findAngleBetweenAngles(start_angle_in_radians, end_angle_in_radians, direction, angles_to_test_in_radians, (fig_num))
%
% INPUTS:
%
%      start_angle_in_radians: the start angle in radians
%
%      end_angle_in_radians: the end angle in radians
%
%      direction: the direction connecting the start and end angles to
%      check. Enter 1 for clockwise, anything else for counter-clockwise.
%
%      angles_to_test: a vector of [N x 1] angles in radians to check
%
%      (OPTIONAL INPUTS)
%
%      fig_num: a figure number to plot results.
%
% OUTPUTS:
%
%      isAngleBetween: an [Nx1] vector containing the 1 if the respective test angle
%      is between the start and end angles, 0 if not
%
% DEPENDENCIES:
%
%      fcn_DebugTools_checkInputsToFunctions
%      fcn_geometry_circleCenterFrom3Points
%      fcn_geometry_arcDirectionFrom3Points
%      fcn_geometry_findAngleUsing2PointsOnCircle
%
% EXAMPLES:   
%
% See the script: script_test_fcn_geometry_findAngleBetweenAngles
% for a full test suite.
%
% This function was written on 2024_01_07 by S. Brennan
% Questions or comments? sbrennan@psu.edu 

% Revision history:
% 2024_01_07 - sbrennan@psu.edu
% -- original write of the code


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
        narginchk(4,5);

        % Check the start_angle_in_radians input
        fcn_DebugTools_checkInputsToFunctions(...
            start_angle_in_radians, '1column_of_numbers',[1 1]);

        % Check the end_angle_in_radians input
        fcn_DebugTools_checkInputsToFunctions(...
            end_angle_in_radians, '1column_of_numbers',[1 1]);

        % Check the direction input
        fcn_DebugTools_checkInputsToFunctions(...
            direction, '1column_of_numbers',[1 1]);

        % angles_to_test_in_radians the direction input
        fcn_DebugTools_checkInputsToFunctions(...
            angles_to_test_in_radians, '1column_of_numbers',[1 2]);
    end
end

% Does user want to show the plots?
flag_do_plots = 0;
if (5 == nargin) && (0==flag_max_speed)
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

% Subtract the start angle from all other angles, effectively rotating the
% origin to the start angle (simplifies calculations)
shifted_end_angle_in_radians = end_angle_in_radians - start_angle_in_radians;
shifted_angles_to_test_in_radians = angles_to_test_in_radians - start_angle_in_radians;

% Perform modulo operation to convert all angles to 0 to 2pi range
modulo_shifted_end_angle_in_radians = mod(shifted_end_angle_in_radians,2*pi);
modulo_shifted_angles_to_test_in_radians = mod(shifted_angles_to_test_in_radians,2*pi);

% Check conditions based on direction
if direction==1
    isAngleBetween = modulo_shifted_angles_to_test_in_radians<=modulo_shifted_end_angle_in_radians + eps*100;
else
    isAngleBetween = modulo_shifted_angles_to_test_in_radians>=modulo_shifted_end_angle_in_radians - eps*100;
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
if flag_do_plots

    temp_h = figure(fig_num);
    flag_rescale_axis = 0;
    if isempty(get(temp_h,'Children'))
        flag_rescale_axis = 1;
    end        

    if flag_rescale_axis == 1
        axis equal;
        hold on;
        grid on;

        % Plot the unit circle
        fcn_geometry_plotCircle([0 0],1);

        % Plot the input angles
        current_figure.h1 = plot([0; cos(start_angle_in_radians)], [0; sin(start_angle_in_radians)],'g.-','MarkerSize',20,'LineWidth',3);
        current_figure.h2 = plot([0; cos(end_angle_in_radians)], [0; sin(end_angle_in_radians)],'r.-','MarkerSize',20,'LineWidth',3);

        % Plot the angles_to_test_in_radians
        current_figure.h3 = plot(cos(angles_to_test_in_radians),sin(angles_to_test_in_radians),'k.','MarkerSize',20,'LineWidth',3);

        % Plot the results of isAngleBetween
        current_figure.h4 = plot(cos(angles_to_test_in_radians(isAngleBetween)),sin(angles_to_test_in_radians(isAngleBetween)),'c.','MarkerSize',10,'LineWidth',3);

        % Make axis slightly larger?
        if flag_rescale_axis
            temp = axis;
            %     temp = [min(points(:,1)) max(points(:,1)) min(points(:,2)) max(points(:,2))];
            axis_range_x = temp(2)-temp(1);
            axis_range_y = temp(4)-temp(3);
            percent_larger = 0.3;
            axis([temp(1)-percent_larger*axis_range_x, temp(2)+percent_larger*axis_range_x,  temp(3)-percent_larger*axis_range_y, temp(4)+percent_larger*axis_range_y]);
        end


        % Save results
        set(gcf,'UserData',current_figure);
    else % Update the results

        % Grab plot information
        current_figure = get(fig_num,'UserData');

        % Update the input angles
        set(current_figure.h1,'XData', [0; cos(start_angle_in_radians)]);
        set(current_figure.h1,'YData', [0; sin(start_angle_in_radians)]);
        set(current_figure.h2,'XData', [0; cos(end_angle_in_radians)]);
        set(current_figure.h2,'YData', [0; sin(end_angle_in_radians)]);
        set(current_figure.h3,'XData', cos(angles_to_test_in_radians));
        set(current_figure.h3,'YData', sin(angles_to_test_in_radians));
        set(current_figure.h4,'XData', cos(angles_to_test_in_radians(isAngleBetween)));
        set(current_figure.h4,'YData', sin(angles_to_test_in_radians(isAngleBetween)));

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





