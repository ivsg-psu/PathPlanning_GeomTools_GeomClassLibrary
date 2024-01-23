function fcn_geometry_plotSphere(centers,radii,varargin)   
% fcn_geometry_plotSphere -  plots a sphere by creating a vector of angles
% spaced 0.01 radians apart, and plotting this as a line around the
% perimeter.
%
% FORMAT:
%
%     fcn_geometry_plotSphere(...
%     centers,...
%     radii,...
%     (fig_num))
%
% INPUTS:
%
%      centers: an [N x 2] vector in [x y] of the points of sphere centers
%
%      radii: a [N x 1] vector of the radii of the spheres (to avoid
%      calculation time)
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
%      fcn_DebugTools_checkInputsToFunctions
%
% EXAMPLES:
%
% See the script: script_test_fcn_geometry_plotSphere
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
    narginchk(2,4); 
    
    % Check the centers input
    fcn_DebugTools_checkInputsToFunctions(...
        centers, '3column_of_numbers');
    
    % Use number of radii to calculate the number of centers
    Nspheres = length(centers(:,1));
    
    % Check the radii input
    fcn_DebugTools_checkInputsToFunctions(...
        radii, '1column_of_numbers',Nspheres);
    
end
    

% Does user want to specify the figure?
flag_do_plot = 1;
fig_num = [];
if 4 == nargin
    fig_num = varargin{end};
    flag_do_plot = 1;
else
    if flag_do_debug
        fig = figure;
        fig_num = fig.Number;
        flag_do_plot = 1;
        flag_new_figure = 1;
    end
end
if isempty(fig_num)
    fig = gcf;
    fig_num = fig.Number;
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



% Set plotting defaults
plot_str = 'b-';
plot_type = 1;  % Plot type refers to 1: a string is given or 2: a color is given - default is 1

% Check to see if user passed in a string or color style?
if 3 <= nargin
    input = varargin{1};
    if ~isempty(input)
        plot_str = input;
        if isnumeric(plot_str)  % Numbers are a color style
            plot_type = 2;
        end
    end
end

% Get a default sphere construct
[X,Y,Z] = sphere;


% Fix the size and position
X2 = X * radii + centers(1);
Y2 = Y * radii + centers(2);
Z2 = Z * radii + centers(3);

% Plot the results in point space
temp_h = figure(fig_num);
flag_rescale_axis = 0;
if isempty(get(temp_h,'Children'))
    flag_rescale_axis = 1;
    hold on;
    axis equal
    grid minor;
else
    hold on;
    axis equal
end

% % Get the color ordering?
% try
%     color_ordering = orderedcolors('gem12');
% catch
%     color_ordering = colororder;
% end
% 
% N_colors = length(color_ordering(:,1));


% Make plots
if plot_type==1
    h_plot = surf(X2,Y2,Z2,'EdgeColor',[0 0.4 0],'FaceAlpha',0.1,'EdgeAlpha',0.1);
    % colormap([0 1 0]);
    %  plot(x_sphere,y_sphere,plot_str);
elseif plot_type==2
    h_plot = surf(X2,Y2,Z2);
    %  plot(x_sphere,y_sphere,'Color',plot_str);
end

view(3);

% Make axis slightly larger?
if flag_rescale_axis
    temp = axis;
    %     temp = [min(points(:,1)) max(points(:,1)) min(points(:,2)) max(points(:,2))];
    axis_range_x = temp(2)-temp(1);
    axis_range_y = temp(4)-temp(3);
    percent_larger = 0.3;
    axis([temp(1)-percent_larger*axis_range_x, temp(2)+percent_larger*axis_range_x,  temp(3)-percent_larger*axis_range_y, temp(4)+percent_larger*axis_range_y]);
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
