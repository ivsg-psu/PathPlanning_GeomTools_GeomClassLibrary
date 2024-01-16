function arc_points = fcn_geometry_plotArc(centers, radii, start_angle_in_radians, end_angle_in_radians, varargin)   
% fcn_geometry_plotArc -  plots an arc by creating a vector of angles
% spaced a fixed angle apart, and plotting this from the start angle to the
% end angle. 
% 
% NOTE: If the end angle is numerically larger than the start angle, the
% plot will be in the positive direction; otherwise, the plot will be in
% the negative direction starting at the start angle, proceeding to the end
% angle.
% 
% NOTE: The user can specify the fixed angle used to space the plotting
% points via an optional input, degree_step.
%
% FORMAT:
%
%     arc_points = fcn_geometry_plotArc(...
%     centers,...
%     radii,...
%     start_angle_in_radians, 
%     end_angle_in_radians,
%     (degree_step),
%     (format),
%     (fig_num))
%
% INPUTS:
%
%      centers: an [N x 2] vector in [x y] of the points of circle centers
%
%      radii: a [N x 1] vector of the radii of the circles (to avoid
%      calculation time)
%
%      start_angle_in_radians: the starting angle of the arc, in radians
% 
%      end_angle_in_radians: the starting angle of the arc, in radians
%
%      (OPTIONAL INPUTS)
%
%      degree_step: how many degrees between plotting points. Default is 1
%      degree.
%
%      format:
%
%        A format string, e.g. 'b-', that dictates the plot style or
%        A color vector, e.g. [1 0 0.23], that dictates the line color
%
%      fig_num: a figure number to plot results.
%
% OUTPUTS:
%
%      arc_points: the [x y] coordinates of the arc points. If N
%      centers and radii are given, with N>1, then arc_points will be a
%      cell array of points.
%
% DEPENDENCIES:
%
%      fcn_geometry_checkInputsToFunctions
%
% EXAMPLES:
%
% See the script: script_test_fcn_geometry_plotArc
% for a full test suite.
%
% This function was written on 2020_10_13 by S. Brennan
% Questions or comments? sbrennan@psu.edu

% Revision History:
% 2021-05-22
% -- new function from fcn_geometry_findAngleUsing3PointsOnCircle
% -- eliminates repo on fcn_plotCircles
% 2024_01_16
% -- updated comments
% -- added fast mode by allowing fig_num set to -1
% -- added degree_step as an optional input
% -- fixed bug in figure argument input
% -- added arc_points as outputs


%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==7 && isequal(varargin{end},-1))
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

%% check input arguments?
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

if flag_check_inputs    
    % Are there the right number of inputs?
    narginchk(4,7); 
    
    % Check the centers input
    fcn_geometry_checkInputsToFunctions(...
        centers, '2column_of_numbers');
    
    % Use number of radii to calculate the number of centers
    Ncircles = length(centers(:,1));
    
    % Check the radii input
    fcn_geometry_checkInputsToFunctions(...
        radii, '1column_of_numbers',Ncircles);
    
end

% Does user want to specify the degree_step?
degree_step = 1;
if 5 == nargin
    temp = varargin{1};
    if ~isempty(temp)
        degree_step = 1;
    end
end

% Set plotting defaults
plot_str = 'b-';
plot_type = 1;  % Plot type refers to 1: a string is given or 2: a color is given - default is 1

% Check to see if user passed in a string or color style?
if 6 <= nargin
    input = varargin{2};
    if ~isempty(input)
        plot_str = input;
        if isnumeric(plot_str)  % Numbers are a color style
            plot_type = 2;
        end
    end
end

% Does user want to specify the figure?
flag_do_plot = 1;
if (0==flag_max_speed) && (7==nargin)
    temp = varargin{end};
    if ~isempty(temp)
        fig_num = temp;
        flag_do_plot = 1;
    else
        flag_do_plot = 0;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _
%  |  \/  |     (_)
%  | \  / | __ _ _ _ __
%  | |\/| |/ _` | | '_ \
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Set angles for plotting
if start_angle_in_radians>end_angle_in_radians
    degree_step = -1*degree_step;
end

% Loop through the arcs, prepping data for plotting each
x_arc = [];
y_arc = [];
for ith_arc = 1:Ncircles 
    angles = (start_angle_in_radians(ith_arc,1):degree_step*pi/180:end_angle_in_radians(ith_arc,1))';
    % Nangles = length(angles(:,1));

    xdata = centers(ith_arc,1)+radii(ith_arc)*cos(angles);
    ydata = centers(ith_arc,2)+radii(ith_arc)*sin(angles);

    x_arc = [x_arc; NaN; xdata]; %#ok<AGROW>
    y_arc = [y_arc; NaN; ydata]; %#ok<AGROW>

    if Ncircles==1
        arc_points = [x_arc y_arc];
    else
        arc_points{ith_arc} = [x_arc y_arc];
    end
end



%% Plot results?
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
    axis equal
    grid minor;
    xlabel('X [meters]');
    ylabel('Y [meters]')

    for ith_arc = 1:Ncircles
        if Ncircles==1
            x_arc = arc_points(:,1);
            y_arc = arc_points(:,2);
        else
            x_arc = arc_points{ith_arc}(:,1);
            y_arc = arc_points{ith_arc}(:,2);
        end

        % Make plots
        if plot_type==1
            if length(plot_str)<=2
                plot(x_arc,y_arc,plot_str);
            else
                plot(x_arc,y_arc,eval(plot_str));
            end
        elseif plot_type==2
            plot(x_arc,y_arc,'Color',plot_str);
        end
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
end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file); 
end
