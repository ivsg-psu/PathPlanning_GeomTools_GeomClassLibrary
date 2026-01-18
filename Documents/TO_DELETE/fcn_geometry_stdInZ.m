function [mean_std_in_z_driven_path,mean_std_in_z_not_driven_path,max_std_in_z_not_driven_path]...
    = fcn_geometry_stdInZ(LiDAR_allPoints,updated_original_mapped_grids,gridIndices,gridIndices_cell_array,...
    gridCenters_updated_original_mapped_grids, updated_current_mapped_grids, current_grid_numbers_of_driven_path, ...
    gridCenters_driven_path,grid_AABBs, grid_size,varargin)

% This is the functionalized part of Step 8: Standard deviation in Z
%
% FORMAT: 
% function [] = fcn_geometry_stdInZ(LiDAR_allPoints,original_mapped_grids,grid_Indices_cell_array,...
%    gridCenters_mapped_grids,gridCenters_driven_path,grid_AABBs, grid_size,varargin)
%
% INPUTS:
% LiDAR_allPoints
% updated_original_mapped_grids
% gridIndices
% gridIndices_cell_array
% gridCenters_updated_original_mapped_grids
% updated_current_mapped_grids
% current_grid_numbers_of_driven_path
% gridCenters_driven_path
% grid_AABBs
% grid_size
%
% OPTIONAL INPUTS:
% fig_num_ENU_statistic_three
% fig_num_1
% fig_num_2
%
% OUTPUTS:
%
% DEPENDENCIES:
%
% EXAMPLES:
% See the script: fcn_geometry_stdInZ
% 
% Code written by:
% Aneesh Batchu 
% 
% Revision History -- Aleksandr Goncharov
% 7/29/2024 - Functionalized the Code
%
%% Debug and Max speed
% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.

flag_max_speed = 0;
if (nargin==13 && isequal(varargin{end},-1))
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
        narginchk(10,13);
 
    end
end

%Figure fig_num_ENU_statistic_three
%ENU_3D_fig_num 
fig_num_ENU_statistic_three = [];
if (11<=nargin)
    temp = varargin{1};
    if ~isempty(temp)
        fig_num_ENU_statistic_three = temp;
    end
end

%Figure fig_num_1
%Figure fig_num_ENU_statistic_three
%ENU_3D_fig_num 
fig_num_1 = [];
if (12<=nargin)
    temp = varargin{2};
    if ~isempty(temp)
        fig_num_1 = temp;
    end
end

%Figure fig_num_2
fig_num_2 = []; % Default is to have no figure
flag_do_plots = 0;
if (0==flag_max_speed) && (13<= nargin)
    temp = varargin{end};
    if ~isempty(temp)
        fig_num_2 = temp;
        flag_do_plots = 1;
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

input_points =  LiDAR_allPoints(:,1:3);

% The indices of the mapped grids are extracted and concatenated
original_mapped_gridIndices_cell = gridIndices_cell_array(updated_original_mapped_grids);

% Total number of mapped grids
total_mapped_grids = length(original_mapped_gridIndices_cell);

% Standard deviations in orthogonal distances of points in the grid to
% plane
% standard_deviation_in_plane_orthogonals = zeros(total_mapped_grids,1);
standard_deviation_in_z = zeros(total_mapped_grids,1);

for ith_mapped_grid = 1:total_mapped_grids
    % points = input_points(original_mapped_gridIndices_cell{ith_mapped_grid},:);
    % points = points(~isnan(points(:,1)),:);
    [~, standard_deviation_in_z(ith_mapped_grid,:), ~, ~, ~, ~] =...
        fcn_geometry_fitPlaneLinearRegression(input_points(original_mapped_gridIndices_cell{ith_mapped_grid},:),-1);
    % mean_z_of_mapped_grids(ith_mapped_grid,:) = mean(input_points(original_mapped_gridIndices_cell{ith_mapped_grid},3));
    % z_diff_mapped_grids = abs(diff(mean_z_of_mapped_grids));
end

figure(fig_num_ENU_statistic_three);clf;

hold on
axis on
grid on
xlabel('X[m]')
ylabel('Y[m]')
title('Grid centers mapped grids in ENU')

% % plot(gridCenters_required_point_density(:,1), gridCenters_required_point_density(:,2), '.','MarkerSize',40,'Color',[0.2 0.2 0.2]);
% plot(gridCenters_drivable_grids(:,1), gridCenters_drivable_grids(:,2), '.','MarkerSize',50,'Color',[0 1 0]);
% plot(gridCenters_non_drivable_grids(:,1), gridCenters_non_drivable_grids(:,2), '.','MarkerSize',50,'Color',[1 0 0]);

plot(gridCenters_updated_original_mapped_grids(:,1), gridCenters_updated_original_mapped_grids(:,2), '.','MarkerSize',45,'Color',[0.2 0.2 0.2]);

% plot the grids in the driven path
plot(gridCenters_driven_path(:,1),gridCenters_driven_path(:,2),'o','MarkerSize',10,'Color',[0 1 0], 'LineWidth',2) % points strictly inside


for ith_text = 1:length(updated_current_mapped_grids(:,1))
    current_text = sprintf('%.0d',updated_current_mapped_grids(ith_text));
    % Place the text on the grid center
    text(gridCenters_updated_original_mapped_grids(ith_text,1), gridCenters_updated_original_mapped_grids(ith_text,2),current_text,'Color',[1 1 1],'HorizontalAlignment','center','FontSize', 6, 'FontWeight','bold');
end


if flag_do_plots
    % Plot grid lines and standard deviation
    figure(fig_num_1)

    hold on
    grid on
    xlabel('X[m]')
    ylabel('Y[m]')
    title('Standard deviation of Z of mapped grids')
end
% Pre-allocation: To find the total number of points in each grid
total_points_in_mapped_grids = zeros(length(updated_original_mapped_grids),1);

% Pre-allocation: these are the coordinates of the corners of each grid
% line. Total no of lines required for each grid including NaNs is 11. Therefore, 11 is multiplied
gridlines_mapped_grids = zeros(11*length(updated_original_mapped_grids),2); % length(gridlines) = 11

for ith_domain = 1:length(updated_original_mapped_grids)
    % Get current color
    % current_color = fcn_geometry_fillColorFromNumberOrName(ith_domain);

    current_color = [0.5 0.5 0.5];

    if flag_do_plots
        % Plot current AABB
        current_AABB = grid_AABBs(updated_original_mapped_grids(ith_domain),1:4);
    end

    % Nudge the current AABB inward
    current_AABB = current_AABB + grid_size/100*[1 -1 1 -1];

    % Calculate the gridlines
    gridlines = [...
        current_AABB(1,1) current_AABB(1,3); ...
        current_AABB(1,1) current_AABB(1,4); ...
        nan nan;
        current_AABB(1,2) current_AABB(1,3); ...
        current_AABB(1,2) current_AABB(1,4); ...
        nan nan;
        current_AABB(1,1) current_AABB(1,3); ...
        current_AABB(1,2) current_AABB(1,3); ...
        nan nan;
        current_AABB(1,1) current_AABB(1,4); ...
        current_AABB(1,2) current_AABB(1,4); ...
        ];


    % Get all points in this domain and plot them
    rows_in_domain = gridIndices == updated_original_mapped_grids(ith_domain);
    points_in_domain = input_points(rows_in_domain,:);

    total_points_in_mapped_grids(ith_domain) = length(points_in_domain);

    length_gridlines = length(gridlines);
    % % Plot the mapped points green
    % p1 = plot(points_in_original_mapped_grids(:,1),points_in_original_mapped_grids(:,2),'.','MarkerSize',20,'Color',[0.4660 0.6740 0.1880]);

    if flag_do_plots
        % Plot the result
        % p1 = plot(gridlines(:,1),gridlines(:,2),'-','Color',current_color,'LineWidth',3);
        plot(gridlines(:,1),gridlines(:,2),'-','Color',current_color,'LineWidth',3);
        % hold on
        % plot(points_in_domain(:,1),points_in_domain(:,2),'.','Color',[0.2,0.2,0.2],'MarkerSize',10)
        gridlines_mapped_grids(1+(ith_domain-1)*length_gridlines:ith_domain*length_gridlines,:) = gridlines;
    end

end

% plot the grids in the driven path
plot(gridCenters_driven_path(:,1),gridCenters_driven_path(:,2),'o','MarkerSize',30,'Color',[0 1 0], 'LineWidth',2) % points strictly inside

% [0.4660 0.6740 0.1880] - Green (Not so bright)
standard_deviation_in_z_round = round(standard_deviation_in_z,3);
% mapped_grid_numbers = 1:length(gridCenters_required_point_density(:,1));
for ith_text = 1:length(gridCenters_updated_original_mapped_grids(:,1))
    current_text = sprintf('%.3f',standard_deviation_in_z_round(ith_text));
    % Place the text on the grid center
    text(gridCenters_updated_original_mapped_grids(ith_text,1), gridCenters_updated_original_mapped_grids(ith_text,2),current_text,'Color',[0 0 0],'HorizontalAlignment','center','FontSize', 10, 'FontWeight','bold','FontSmoothing','on');
end


% Find driven path indices in current mapped grids
driven_path_grid_indices_in_current_mapped_grids = ismember(updated_current_mapped_grids,current_grid_numbers_of_driven_path);

% Standard deviation in Z of driven path grids
std_in_z_driven_path = standard_deviation_in_z(driven_path_grid_indices_in_current_mapped_grids);

% Standard deviation in Z of other mapped grids
std_in_z_other_mapped_grids = standard_deviation_in_z(~driven_path_grid_indices_in_current_mapped_grids);




% Find mean std in z of driven path
mean_std_in_z_driven_path = mean(std_in_z_driven_path);

% Find mean std in z of not driven path
mean_std_in_z_not_driven_path = mean(std_in_z_other_mapped_grids(~isnan(std_in_z_other_mapped_grids)));

% Find max std in z of not driven path
max_std_in_z_not_driven_path = max(std_in_z_other_mapped_grids);

% Std Threshold
% Instead of choosing 6, try to find a ratio
% ratio: mean_std_in_z_driven_path/mean_std_of_all_grids
% std_threshold = mean_std_in_z_driven_path*6;
%
% disp(std_threshold)

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


%Plotting the Results

if flag_do_plots


    % Plot grid lines and standard deviation
    figure(fig_num_2);
    clf;

    hold on
    grid on
    xlabel('Mapped grid centers')
    ylabel('Standard deviation in Z')
    title('Mapped grid centers vs standard deviation in Z ')

    plot(updated_current_mapped_grids, standard_deviation_in_z,'.','MarkerSize',30,'Color',[0.2 0.2 0.2])
    plot(current_grid_numbers_of_driven_path, std_in_z_driven_path,'o','MarkerSize',10,'Color',[0 1 0], 'LineWidth',1.5)
    plot(updated_current_mapped_grids(~driven_path_grid_indices_in_current_mapped_grids), std_in_z_other_mapped_grids,'.','MarkerSize',10,'Color',[1 0 0])



end
end