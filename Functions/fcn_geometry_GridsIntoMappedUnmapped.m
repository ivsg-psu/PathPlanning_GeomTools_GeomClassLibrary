function [original_grids_with_required_point_density,gridCenters_low_point_density,gridCenters_required_point_density,current_grids_with_low_point_density,current_grids_with_required_point_density] = fcn_geometry_GridsIntoMappedUnmapped(point_density,total_N_points_in_each_grid,grids_greater_than_zero_points,gridCenters,varargin)
%% fcn_geometry_GridsIntoMappedUnmapped 
% classify grids with more than zero points into mapped and unmapped  
% 
% FORMAT:
%
%      [original_grids_with_required_point_density] = fcn_geometry_classifyGridsintoMappedUnmapped(point_density,total_N_points_in_each_grid,grids_greater_than_zero_points,grid_centers,varargin)
%
% INPUTS:     
%       
%      point_density: density of point
%
%      total_N_points_in_each_grid
%
%      grids_greater_than_zero_points: grids with more than zero points
%
%      gridCenters: grid Centers of the grids 
%
%      (OPTIONAL INPUTS)
%
%      fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed.
%
% OUTPUTS:
%
%      original_grids_with_required_point_density: The grid numbers belong to the initial grid indices cell array
%
%      gridCenters_low_point_density: Grid Centers of the grids with low point density (Unmapped grid centers)
%
%      gridCenters_required_point_density: Grid Centers of the grids with required point density (Mapped grid centers)
%
%      current_grids_with_low_point_density: Current grid numbers of the grids with low point density
%
%      current_grids_with_required_point_density: grid numbers of the grids with respect to grids_greater_than_zero_points
%
% DEPENDENCIES:
%
%      (none)
%
% EXAMPLES:
%
%       See the script:
%       script_test_fcn_geometry_classifyGridsintoMappedUnmapped
%       test suite.
%
% Revision History
% 2024_06_15 - Aneesh Batchu
% -- Wrote the code originally
% 2024_07_15 - Aneesh Batchu
% -- Seperated this code from fcn_geometry_surfaceAnalysis
% 2024_07_15 - Jiabao Zhao
% -- Functionalize the code

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==5 && isequal(varargin{end},-1))
    flag_do_debug = 0; % % % % Flag to plot the results for debugging
    flag_check_inputs = 0; % Flag to perform input checking
    flag_max_speed = 1;
else
    % Check to see if we are externally setting debug mode to be "on"
    flag_do_debug = 0; % % % % Flag to plot the results for debugging
    flag_check_inputs = 1; % Flag to perform input checking
    MATLABFLAG_GEOMETRY_FLAG_CHECK_INPUTS = getenv("MATLABFLAG_GEOMETRY_FLAG_CHECK_INPUTS");
    MATLABFLAG_GEOMETRY_FLAG_DO_DEBUG = getenv("MATLABFLAG_GEOMETRY_FLAG_DO_DEBUG");
    if ~isempty(MATLABFLAG_GEOMETRY_FLAG_CHECK_INPUTS) && ~isempty(MATLABFLAG_GEOMETRY_FLAG_DO_DEBUG)
        flag_do_debug = str2double(MATLABFLAG_GEOMETRY_FLAG_DO_DEBUG); 
        flag_check_inputs  = str2double(MATLABFLAG_GEOMETRY_FLAG_CHECK_INPUTS);
    end
end

% flag_do_debug = 1;

if flag_do_debug
    st = dbstack; %#ok<*UNRCH>
    fprintf(1,'STARTING function: %s, in file: %s\n',st(1).name,st(1).file);
    debug_fig_num = 999978; %#ok<NASGU>
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
        narginchk(1,5);

        % % Check the points input to be length greater than or equal to 2
        % fcn_DebugTools_checkInputsToFunctions(...
        %     points, '2column_of_numbers',[2 3]);
        %
        % % Check the transverse_tolerance input is a positive single number
        % fcn_DebugTools_checkInputsToFunctions(transverse_tolerance, 'positive_1column_of_numbers',1);
        %
        % % Check the station_tolerance input is a positive single number
        % if ~isempty(station_tolerance)
        %     fcn_DebugTools_checkInputsToFunctions(station_tolerance, 'positive_1column_of_numbers',1);
        % end
    end
end

% Does user want to specify fig_num?
flag_do_plots = 0;
if 2<= nargin && 0==flag_max_speed
    temp = varargin{end};
    if ~isempty(temp)
        fig_num = temp;
        flag_do_plots = 1;
    end
end

%% Solve for the Maxs and Mins
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _       
%  |  \/  |     (_)      
%  | \  / | __ _ _ _ __  
%  | |\/| |/ _` | | '_ \ 
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find grids with low point density but not zero point density (Unmapped grid centers)
% Original: The grid numbers belong to the initial grid indices cell array
% The numbering does not start from "1"
original_grids_with_low_point_density = grids_greater_than_zero_points((total_N_points_in_each_grid(grids_greater_than_zero_points,1) > 0) & (total_N_points_in_each_grid(grids_greater_than_zero_points,1) < point_density)); 

% Find mapped grids 
% Mapped grids. Later these grids are classified into drivable and non
% drivable
% Original: The grid numbers belong to the initial grid indices cell array
% The numbering does not start from "1"
original_grids_with_required_point_density = grids_greater_than_zero_points(total_N_points_in_each_grid(grids_greater_than_zero_points,1) >= point_density); 

% Grid Centers of the grids with low point density (Unmapped grid centers)
gridCenters_low_point_density = gridCenters(original_grids_with_low_point_density,1:3);  

% Grid Centers of the grids with required point density (Mapped grid centers)
gridCenters_required_point_density = gridCenters(original_grids_with_required_point_density,1:3);  

% Current grid numbers of the grids with low point density
% These are the grid numbers of the grids with respect to grids_greater_than_zero_points
% The numbering starts with "1"
current_grids_with_low_point_density = find((total_N_points_in_each_grid(grids_greater_than_zero_points,1) > 0) & (total_N_points_in_each_grid(grids_greater_than_zero_points,1) < point_density)); 

% Current grid numbers of the grids with required point density
% These are the grid numbers of the grids with respect to grids_greater_than_zero_points
% The numbering starts with "1"
current_grids_with_required_point_density = find(total_N_points_in_each_grid(grids_greater_than_zero_points,1) >= point_density); 

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
    fig_num = 11; 
    figure(fig_num);clf

    marker_size = 30;
    RGB_triplet = [0.8, 0.8, 0.8];
    legend_option = 1;
    legend_name = 'Unmapped grids';
    legend_position = [];
    [~] = fcn_geometry_plotGridCenters(gridCenters_low_point_density,marker_size,RGB_triplet,legend_option,legend_name,legend_position,fig_num);

    % plot grid centers
    marker_size = 30;
    RGB_triplet = [0.2, 0.2, 0.2];
    legend_option = 1;
    legend_name = 'Mapped grids';
    legend_position = [];
    [~] = fcn_geometry_plotGridCenters(gridCenters_required_point_density,marker_size,RGB_triplet,legend_option,legend_name,legend_position,fig_num);
end % Ends check if plotting

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file);
end
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
