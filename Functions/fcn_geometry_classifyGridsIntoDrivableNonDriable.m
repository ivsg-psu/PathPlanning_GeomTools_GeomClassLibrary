function [standard_deviation_in_z,angle_btw_unit_normals_and_vertical,original_drivable_grids,original_non_drivable_grids,current_drivable_grid_numbers_in_mapped_grids,current_non_drivable_grid_numbers_in_mapped_grids,gridCenters_drivable_grids,gridCenters_non_drivable_grids] = fcn_geometry_classifyGridsIntoDrivableNonDriable(original_grids_with_required_point_density,input_points,std_threshold,theta_threshold,gridCenters, varargin)
%% fcn_geometry_classifyMappedGridsIntoDrivableAndNondrivable
% classify mapped grids into drivable and non-drivable
% 
% FORMAT:
%
%      [standard_deviation_in_z,angle_btw_unit_normals_and_vertical,original_drivable_grids,original_non_drivable_grids,current_drivable_grid_numbers_in_mapped_grids,current_non_drivable_grid_numbers_in_mapped_grids,gridCenters_drivable_grids,gridCenters_non_drivable_grids] = fcn_geometry_classifyMappedGridsIntoDrivableAndNondrivable(original_grids_with_required_point_density,input_points,std_threshold,theta_threshold,gridCenters, varargin)
%
% INPUTS:     
%       
%      original_grids_with_required_point_density:The grid numbers belong to the initial grid indices cell array
%
%      input_points: outer range of mapping data 
%
%      std_threshold: standard deviation of threshold
%
%      theta_threshold: angle of threshold 
%
%      gridCenters: center of grids
%      
%      (OPTIONAL INPUTS)
%
%      fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed.
%
% OUTPUTS:
%
%     standard_deviation_in_z: standard deviation in z direction
% 
%     angle_btw_unit_normals_and_vertical: the angle between unit normals and vertical
% 
%     original_drivable_grids: drivable grids
% 
%     original_non_drivable_grids: non-derivable grids
% 
%     current_drivable_grid_numbers_in_mapped_grids: drivable grid numbers of the mapped grids
% 
%     current_non_drivable_grid_numbers_in_mapped_grids: non-drivable grid numbers of the mapped grids
% 
%     gridCenters_drivable_grids: Grid centers of drivable grids 
% 
%     gridCenters_non_drivable_grids: Grid centers of non-drivable grids 
%
% DEPENDENCIES:
%
%      (none)
%
% EXAMPLES:
%
%       See the script:
%       script_test_fcn_geometry_classifyMappedGridsIntoDrivableAndNondrivable.m for a full
%       test suite.
%
% This function was written on 2024_07_15 by Aneesh Batchu
% -- Seperated this code from fcn_geometry_surfaceAnalysis
% Funclionalize this code on 7/16/2024 by Jiabao Zhao

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==6 && isequal(varargin{end},-1))
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
        narginchk(1,6);

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



% The indices of the mapped grids are extracted and concatenated 
original_mapped_gridIndices_cell = gridIndices_cell_array(original_grids_with_required_point_density); 

% Total number of mapped grids
total_mapped_grids = length(original_mapped_gridIndices_cell); 

% Standard deviations in orthogonal distances of points in the grid to
% plane
% standard_deviation_in_plane_orthogonals = zeros(total_mapped_grids,1); 
standard_deviation_in_z = zeros(total_mapped_grids,1); 

% Unit normal vectors of the plane fits of each mapped grid
unit_normal_vectors = zeros(total_mapped_grids,3); 


% z_height of all the points 
% mean_z_of_mapped_grids = zeros(total_mapped_grids,1); 

% Loop through all the mapped grids, recording standard deviation, unit
% vectors 
% if 0==flag_max_speed
%     h_waitbar = waitbar(0,'Performing surface analysis...');
% end

for ith_mapped_grid = 1:total_mapped_grids
    [~, standard_deviation_in_z(ith_mapped_grid,:), ~, unit_normal_vectors(ith_mapped_grid,:), ~, ~] =...
    fcn_geometry_fitPlaneLinearRegression(input_points(original_mapped_gridIndices_cell{ith_mapped_grid},:),-1);
    mean_z_of_mapped_grids(ith_mapped_grid,:) = mean(input_points(original_mapped_gridIndices_cell{ith_mapped_grid},3));
    z_diff_mapped_grids = abs(diff(mean_z_of_mapped_grids)); 
end

% standard_deviation_in_z = round(standard_deviation_in_z,4); 
if ~isempty(std_threshold) && isempty(theta_threshold)

    % STEP 1: Standard deviation of the orthogonal (perpendicular) distances of
    % the points to the plane (after fit)
    % Find the grids that are within standard deviation limit
    % This is not enough (delta Y) is also important
    % mapped_grids_within_std_threshold = standard_deviation_in_plane_orthogonals < std_threshold;
    mapped_grids_within_std_threshold = standard_deviation_in_z < std_threshold;
    
    % Grids that satisy the conditions of (STEP 1). The grids that
    % are within the std threshold
    mapped_grids_within_all_thresholds = (mapped_grids_within_std_threshold == 1);

    % The angle between unit vertical and the unit_normal_vector is computed to
    % determine how close the normal vector is to vertical direction. In
    % this case, the angle between unit normals and vertical is empty. 
    angle_btw_unit_normals_and_vertical = [];

elseif isempty(std_threshold) && ~isempty(theta_threshold)

    % STEP 2
    % Comparing normal vector with verticle direction
    unit_vector_vertical_direction = [0 0 1];

    % The dot product is computed to find the angle between the vectors
    dot_product = sum(unit_normal_vectors.*unit_vector_vertical_direction,2);

    % The angle between unit vertical and the unit_normal_vector is computed to
    % determine how close the normal vector is to vertical direction.
    angle_btw_unit_normals_and_vertical = acos(dot_product);

    % Find the grids (with a fitted plane) that are within the vertical
    % threshold (change to a different name: vertical threshold)
    mapped_grids_within_vertical_threshold = angle_btw_unit_normals_and_vertical < theta_threshold;

    % Grids that satisy the conditions (STEP 2). The grids that
    % are within the vertical threshold
    mapped_grids_within_all_thresholds = (mapped_grids_within_vertical_threshold == 1);

else

    % STEP 1: Standard deviation of the orthogonal (perpendicular) distances of
    % the points to the plane (after fit)
    % Find the grids that are within standard deviation limit
    % This is not enough (delta Y) is also important
    % mapped_grids_within_std_threshold = standard_deviation_in_plane_orthogonals < std_threshold;
    mapped_grids_within_std_threshold = standard_deviation_in_z < std_threshold;

    % STEP 2
    % Comparing normal vector with verticle direction
    unit_vector_vertical_direction = [0 0 1];

    % The dot product is computed to find the angle between the vectors
    dot_product = sum(unit_normal_vectors.*unit_vector_vertical_direction,2);

    % The angle between unit vertical and the unit_normal_vector is computed to
    % determine how close the normal vector is to vertical direction.
    angle_btw_unit_normals_and_vertical = acos(dot_product);

    % Find the grids (with a fitted plane) that are within the vertical
    % threshold (change to a different name: vertical threshold)
    mapped_grids_within_vertical_threshold = angle_btw_unit_normals_and_vertical < theta_threshold;

    % Grids that satisy the conditions of (STEP 1 & STEP 2). The grids that
    % are within the standar deviation and vertical threshold
    mapped_grids_within_vertical_and_std_thresholds = (mapped_grids_within_vertical_threshold == 1) & (mapped_grids_within_std_threshold == 1);
    
    % mapped grids within all the thresholds 
    mapped_grids_within_all_thresholds = mapped_grids_within_vertical_and_std_thresholds;

end

% Find the drivable grids (original)
original_drivable_grids = original_grids_with_required_point_density(mapped_grids_within_all_thresholds); 

% Find the non-drivable grids (original)
original_non_drivable_grids = original_grids_with_required_point_density(mapped_grids_within_all_thresholds == 0);

% Final drivable grid numbers of the mapped grids
current_drivable_grid_numbers_in_mapped_grids = find(ismember(original_grids_with_required_point_density, original_drivable_grids));

% Final non drivable grid numbers of the mapped grids
current_non_drivable_grid_numbers_in_mapped_grids = find(ismember(original_grids_with_required_point_density, original_non_drivable_grids));

% Grid centers of drivable grids 
gridCenters_drivable_grids = [gridCenters(original_drivable_grids,1), gridCenters(original_drivable_grids,2), ones(length(original_drivable_grids),1)]; 

% Grid centers of nondrivable grids
gridCenters_non_drivable_grids = [gridCenters(original_non_drivable_grids,1), gridCenters(original_non_drivable_grids,2), zeros(length(original_non_drivable_grids),1)]; 

% Concatenate the grid centers of drivable and non-drivable grids (2D)
% gridCenters_mapped_grids = [gridCenters_drivable_grids; gridCenters_non_drivable_grids];

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
        
marker_size = 30;
RGB_triplet = [0, 1, 0]; 
legend_option = 1;
legend_name = 'Drivable grids';
legend_position = [];
plot_gridCenters_drivable_grids = [gridCenters_drivable_grids(:,1:2), zeros(length(gridCenters_drivable_grids(:,1)),1)]; 
[~] = fcn_geometry_plotGridCenters(plot_gridCenters_drivable_grids,marker_size,RGB_triplet,legend_option,legend_name,legend_position,fig_num);

% plot grid centers
marker_size = 30;
RGB_triplet = [1, 0, 0]; 
legend_option = 1;
legend_name = 'Non-drivable grids';
legend_position = [];
plot_gridCenters_non_drivable_grids = [gridCenters_non_drivable_grids(:,1:2), zeros(length(gridCenters_non_drivable_grids(:,1)),1)]; 
[~] = fcn_geometry_plotGridCenters(plot_gridCenters_non_drivable_grids,marker_size,RGB_triplet,legend_option,legend_name,legend_position,fig_num);
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
