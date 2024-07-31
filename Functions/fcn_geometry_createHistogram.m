function [point_density, counts1, counts2, binEdges] = fcn_geometry_createHistogram(total_points_in_each_grid_in_the_driven_path,...
    total_points_in_each_grid_with_points_greater_than_zero,grid_size, varargin)
%% fcn_geometry_createHistogram
% Create a histogram that shows points per grid
% 
% FORMAT:
%
%      [point_density, counts1, counts2, binEdges] = fcn_geometry_createHistogram(total_points_in_each_grid_in_the_driven_path,...
%      total_points_in_each_grid_with_points_greater_than_zero,grid_size, (fig_num))
%
% INPUTS:     
%       
%      total_points_in_each_grid_in_the_driven_path: The total points in
%      each grid of the driven path of the mapping van
%
%    
%
% OUTPUTS:
%       
%      point_density: the density of point in each grid
%
% DEPENDENCIES:
%
%      (none)
%
% EXAMPLES:
%
%       See the script:
%       script_test_fcn_geometry_createHistogram.m for a full
%       test suite.
%
% This code was written on 2024_07_18 by Aneesh Batchu
% This code was functionalized on 7/29/2024 by Jiabao Zhao
% Questions or comments? jpz5469@psu.edu
%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==4 && isequal(varargin{end},-1))
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
        narginchk(1,4);

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

% Calculate the overlapping counts
% % overlapCounts = min(counts2, counts1);

% Find a ratio
 % point_density = sum(binEdges(index_max_counts2:(index_max_counts2+1)))/2; 

% point_density = floor(mean(total_points_in_each_grid_in_the_driven_path) - 7.5*(std(total_points_in_each_grid_in_the_driven_path)));
% point_density = floor(mean(total_points_in_each_grid_in_the_driven_path) - 1.5*(std(total_points_in_each_grid_in_the_driven_path)));

% point_density = floor(mean(total_points_in_each_grid_in_the_driven_path));


% Minimum number of points required 
point_density = floor(20*((grid_size^2)/(0.3^2))); 
disp('Chosen point density')
disp(point_density)
% mean(total_points_in_each_grid_in_the_driven_path)/mean(total_points_in_each_grid_with_points_greater_than_zero);

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

    hold on
    grid on
    xlabel('Points per grid');
    ylabel('Frequency');
    title('Statistic 1: Determining suitable point density');
    % edges = (floor(min(total_points_in_each_grid_with_points_greater_than_zero)/10)*10):10:(ceil(max(total_points_in_each_grid_with_points_greater_than_zero)/10)*10); % Define the bin edges

    % Create the histogram
    % actual_driving_surface_grids_hist = histogram(total_points_in_each_grid_in_the_driven_path,edges,'Visible','on');
    % total_grids_hist = histogram(total_points_in_each_grid_with_points_greater_than_zero,edges,'Visible','on');
    total_grids_greater_than_zero_hist = histogram(total_points_in_each_grid_with_points_greater_than_zero,20,'Visible','on');
    actual_driven_path_grids_hist = histogram(total_points_in_each_grid_in_the_driven_path,10,'Visible','on');
    % Extract the counts for both histograms
    counts1 = total_grids_greater_than_zero_hist.Values;
    %[~,index_max_counts2] = max(counts1);
    counts2 = actual_driven_path_grids_hist.Values;

    binEdges = total_grids_greater_than_zero_hist.BinEdges;
    disp('Chosen point density')
    disp(point_density)
    % mean(total_points_in_each_grid_in_the_driven_path)/mean(total_points_in_each_grid_with_points_greater_than_zero);

    plot(point_density,0, 'k.', 'MarkerSize',20)
    current_text = sprintf('point density = %.2d',point_density);
    disp('Chosen point density')
    disp(point_density)
    text(650, 200,current_text,'Color',[0 0 0],'HorizontalAlignment','center','FontSize', 12, 'FontWeight','bold');
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
