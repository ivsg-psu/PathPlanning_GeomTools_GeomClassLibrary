function [z,mean_z,deviation_from_mean,angles,mean_angle,angle_difference]...
    = fcn_geometry_plotLidarDeviation(cell_start,cell_end,color_map,varargin)
% Takes in cell_start and cell_end as well as additional variables
% to plot a deviation from the mean in terms of height and angles.
%
% FORMAT: 
% [z,mean_z,deviation_from_mean,angles,mean_angle,angle_difference]...
% = fcn_geometry_plotLidarDeviation(cell_start,cell_end,color_map...
% ,(file_name)(fig1),(fig2),(fig3));
%
% INPUTS:
%       cell_start: starting cell
%       cell_end: ending cell
%       color_map: color map to plot the points in
%       (OPTIONAL INPUTS)
%
%       file_name: if an external file is needed for LiDAR_ENU_cell
%       fig_1: Label for figure 1
%       fig_2: Label for figure 2
%       fig_3: Label for figure 3
%
% OUTPUTS:
%   z, mean_z, deviation_from_mean
%   angles, mean_angle,angle_difference
%
% DEPENDENCIES:
%    The script must either have existing LiDAR_ENU_cell or import a file
%
% EXAMPLES:
%   See the script:
%   script_test_fcn_geometry_plotLidarDeviation
%
%
% Revision History
% Aneesh Batchu - 2024_07_08
% -- wrote the code originally
% Aleksandr Goncharov 2024_07_08
% -- functionalized the code
% Aleksandr Goncharov 2024_07_17
% -- updated the function to the correct format, created a script


%Example:[z,mean_z,deviation_from_mean,angles,mean_angle,angle_difference]...
% = fcn_geometry_plotLidarDeviation(50,59,'jet','TestTrack_Entire_Loop',...
% 111,222,333);

%% Debug and Max speed
% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.

flag_max_speed = 0;
if (nargin==7 && isequal(varargin{end},-1))
    flag_do_debug = 0; % Flag to plot the results for debugging
    flag_max_speed = 1;
else
    % Check to see if we are externally setting debug mode to be "on"
    flag_do_debug = 0; % Flag to plot the results for debugging
    MATLABFLAG_LOADWZ_FLAG_CHECK_INPUTS = getenv("MATLABFLAG_LOADWZ_FLAG_CHECK_INPUTS");
    MATLABFLAG_LOADWZ_FLAG_DO_DEBUG = getenv("MATLABFLAG_LOADWZ_FLAG_DO_DEBUG");
    if ~isempty(MATLABFLAG_LOADWZ_FLAG_CHECK_INPUTS) && ~isempty(MATLABFLAG_LOADWZ_FLAG_DO_DEBUG)
        flag_do_debug = str2double(MATLABFLAG_LOADWZ_FLAG_DO_DEBUG);
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

if flag_max_speed == 0
    % Are there the right number of inputs?
    narginchk(3,7);
end

%file name if needed
if 4<= nargin
    temp = varargin{1};
    if ~isempty(temp)
        file_name=temp;
    end
end

%labeling figures

flag_do_plots = 0;
fig_1=[];
fig_2=[];
fig_3=[];

%figure 1
if flag_max_speed==0
if 5<= nargin
    temp = varargin{2};
    if ~isempty(temp)
        fig_1=temp;
        flag_do_plots=1;
    end
end

%figure 2
if 6<= nargin
    temp = varargin{3};
    if ~isempty(temp)
        fig_2=temp;
        flag_do_plots=1;
    end
end
%figure 3
if 7 == nargin
    temp = varargin{end};
    if ~isempty(temp)
        fig_3=temp;
        flag_do_plots=1;
    end
end
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

LiDAR_data=load(file_name);
LiDAR_ENU_cell = LiDAR_data.LiDAR_ENU_cell;

% LiDAR data
LiDAR_outer_edge = vertcat(LiDAR_ENU_cell{cell_start:cell_end});

x = LiDAR_outer_edge(:,1);
y = LiDAR_outer_edge(:,2);
z = LiDAR_outer_edge(:,3); 

% Calculate the deviation of z-values from the mean
mean_z = mean(z);
deviation_from_mean = z - mean_z;

% Calculate the angles relative to the x-axis
angles = atan2(y, x);

% Calculate the difference between each angle and the mean angle
mean_angle = mean(angles);
angle_difference = angles - mean_angle;

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

if flag_do_plots
    if ~isempty(fig_1)
    figure(fig_1)
    % Create the scatter plot
    scatter3(x, y, z, 20, deviation_from_mean, 'filled');
    colormap(color_map);
    colorbar;
    view(2);
    title('Deviation of z-value from the mean');
    xlabel('x');
    ylabel('y');
    end
    
    if ~isempty(fig_2)
    figure(fig_2)
    % Create the scatter plot
    scatter3(x, y, z, 20, angle_difference, 'filled');
    colormap(color_map);
    colorbar;
    view(2);
    title('Angle difference from the mean angle');
    xlabel('x');
    ylabel('y');
    end
    
    if ~isempty(fig_3)
    figure(fig_3)
    colors = z;
    scatter3(x,y,z,20,colors);
    colormap(color_map); % You can use 'viridis' if you have it, or other colormaps
    colorbar; % Add a color bar to show the color mapping
    view(2)
    end
    
end
end