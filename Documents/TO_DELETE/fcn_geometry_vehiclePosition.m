function [vehicle_change_in_pose_XY,unit_vehicle_change_in_pose_XY,unit_ortho_vehicle_vectors_XY]=...
    fcn_geometry_vehiclePosition(VehiclePose,varargin)

%% fcn_geometry_vehicleOrientation
% Takes in the inputs of where on the track we want to analyze and provides
% a concatenated list of points for that area. 
%
% INPUTS:
% VehiclePose
%
% OPTIONAL INPUTS:
% (ColorTripletTrack)
% (ColorTripletStart)
% (ColorTripletEnd)
% (MarkerSizeTrack)
% (MarkerSizeEnds)
% (fig_num)
%
% REQUIREMENTS: 
% (None)
%
% EXAMPLES:
% script_test_fcn_geometry_vehicleOrientation
%
% 2024_07_17 - S Brennan
% -- wrote the code
% 2024_07_26 - A Goncharov
% -- Functionalized this part of the code
%
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
    if flag_check_inputs
        % Are there the right number of inputs?
        narginchk(1,7);
 
    end
end

%color triplets
ColorTripletTrack=[1 1 0];
if (2<=nargin)
    temp = varargin{1};
    if ~isempty(temp)
        ColorTripletTrack = temp;
    end
end

ColorTripletStart=[0 1 0];
if (3<=nargin)
    temp = varargin{2};
    if ~isempty(temp)
        ColorTripletTrack = temp;
    end
end

ColorTripletEnd= [1 0 0];
if (4<=nargin)
    temp = varargin{3};
    if ~isempty(temp)
        ColorTripletTrack = temp;
    end
end

MarkerSizeTrack = 10;
if (4<=nargin)
    temp = varargin{4};
    if ~isempty(temp)
        MarkerSizeTrack = temp;
    end
end

MarkerSizeEnds = 10;
if (4<=nargin)
    temp = varargin{5};
    if ~isempty(temp)
        MarkerSizeEnds = temp;
    end
end

% Does user want to specify fig_num?
fig_num = []; % Default is to have no figure
flag_do_plots = 0;
if (0==flag_max_speed) && (7<= nargin)
    temp = varargin{end};
    if ~isempty(temp)
        fig_num = temp;
        flag_do_plots = 1;
    end
end


%% Find the grids
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _
%  |  \/  |     (_)
%  | \  / | __ _ _ _ __
%  | |\/| |/ _` | | '_ \
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Calculate the vehicle orientation
vehicle_change_in_pose_XY = diff(VehiclePose(:,1:2));

% Repeat the last value again, since diff removes one row. We want the same
% number of vectors as the number of points, and diff removed one point.
vehicle_change_in_pose_XY = [vehicle_change_in_pose_XY; vehicle_change_in_pose_XY(end,:)];

% Convert these to unit vectors
unit_vehicle_change_in_pose_XY = fcn_geometry_calcUnitVector(vehicle_change_in_pose_XY);

% Find orthogonal vectors by rotating by 90 degrees in the CCW direction
unit_ortho_vehicle_vectors_XY = unit_vehicle_change_in_pose_XY*[0 1; -1 0];



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
    %% Plot the vehicle pose in ENU
    figure(fig_num);

    hold on;
    grid on;
    axis equal

    plot(VehiclePose(:,1),VehiclePose(:,2),'-','Color',[0 0 0],'MarkerSize',30,'LineWidth',3);
    plot(VehiclePose(:,1),VehiclePose(:,2),'-','Color',ColorTripletTrack,'MarkerSize',MarkerSizeTrack,'LineWidth',1);

    % Plot start and end points
    plot(VehiclePose(1,1),VehiclePose(1,2),'.','Color',ColorTripletStart,'MarkerSize',MarkerSizeEnds);
    plot(VehiclePose(end,1),VehiclePose(end,2),'o','Color',ColorTripletEnd,'MarkerSize',MarkerSizeEnds);

    xlabel('East position [m]');
    ylabel('North position [m]');
end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file);
end

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
end