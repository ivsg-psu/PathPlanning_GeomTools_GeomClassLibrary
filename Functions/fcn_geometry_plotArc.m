function fcn_geometry_plotArc(centers, radii, start_angle_in_radians, end_angle_in_radians, varargin)   
% fcn_geometry_plotArc -  plots an arc by creating a vector of angles
% spaced 0.01 radians apart, and plotting this from the start angle to the
% end angle, in the positive direction
%
% FORMAT:
%
%     fcn_geometry_plotArc(...
%     centers,...
%     radii,...
%     start_angle_in_radians, end_angle_in_radians,
%     (fig_num))
%
% INPUTS:
%
%      centers: an [N x 2] vector in [x y] of the points of circle centers
%
%      radii: a [N x 1] vector of the radii of the circles (to avoid
%      calculation time)
%
%      start_angle_in_radians, end_angle_in_radians: the start and end
%      angles respectively of the arc
%
%      (OPTIONAL INPUTS)
%
%      format:
%        A format string, e.g. 'b-', that dictates the plot style or
%        A color vector, e.g. [1 0 0.23], that dictates the line color
%
%      fig_num: a figure number to plot results.
%
% OUTPUTS:
%
%      (none)
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


%% Debugging and Input checks
flag_check_inputs = 1; % Set equal to 1 to check the input arguments
flag_do_plot = 0;      % Set equal to 1 for plotting "extra" figures (not used)
flag_do_debug = 0;     % Set equal to 1 for debugging

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
    narginchk(4,6); 
    
    % Check the centers input
    fcn_geometry_checkInputsToFunctions(...
        centers, '2column_of_numbers');
    
    % Use number of radii to calculate the number of centers
    Ncircles = length(centers(:,1));
    
    % Check the radii input
    fcn_geometry_checkInputsToFunctions(...
        radii, 'column_of_numbers',Ncircles);
    
end
    

% Does user want to specify the figure?
flag_new_figure = 0;
if 6 == nargin
    fig_num = varargin{end};
    flag_do_plot = 1;
else
    if flag_do_debug
        fig = figure;
        fig_num = fig.Number;
        flag_do_plot = 1;
        flag_new_figure = 1;
    else
        fig = gcf;
        fig_num = fig.Number;
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

figure(fig_num);

if flag_new_figure        
    hold on;
    axis equal
    grid minor;
else
    hold on;
    axis equal
end

% Set plotting defaults
plot_str = 'b-';
plot_type = 1;  % Plot type refers to 1: a string is given or 2: a color is given - default is 1

% Check to see if user passed in a string or color style?
if 5 <= nargin
    input = varargin{1};
    if ~isempty(input)
        plot_str = input;
        if isnumeric(plot_str)  % Numbers are a color style
            plot_type = 2;
        end
    end
end


% Set angles for plotting
degree_step = 1;
if start_angle_in_radians>end_angle_in_radians
    degree_step = -1;
end

% Loop through the arcs, prepping data for plotting each
x_circle = [];
y_circle = [];
for ith_arc = 1:Ncircles 
    angles = (start_angle_in_radians(ith_arc,1):degree_step*pi/180:end_angle_in_radians(ith_arc,1))';
    % Nangles = length(angles(:,1));

    xdata = centers(ith_arc,1)+radii(ith_arc)*cos(angles);
    ydata = centers(ith_arc,2)+radii(ith_arc)*sin(angles);

    x_circle = [x_circle; NaN; xdata]; %#ok<AGROW>
    y_circle = [y_circle; NaN; ydata]; %#ok<AGROW>
end

% Make plots
if plot_type==1
    plot(x_circle,y_circle,plot_str);
elseif plot_type==2
    plot(x_circle,y_circle,'Color',plot_str);
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
    % Nothing more to do here!
end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file); 
end
