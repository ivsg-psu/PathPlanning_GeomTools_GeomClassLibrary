function sequence_indicies = fcn_geometry_findPointsInSequence(input_distances, base_point_index, station_tolerance, varargin)
% fcn_geometry_findPointsInSequence
% Finds the indicies before and after a base point that are within a
% tolerance distance sequence of the base point. For example, if the number
% 7 (index 6 below) is the base point and the tolerance is a distance of 2,
% then the inputs of:
%
%     -100 -99 0 3 6 7 8.5 9 22
%
% Returns the indicies: [5 6 7 8] because the numbers from 6 to 9 are all
% in sequence within a tolerance distance spacing of 2. Note: the input
% distances can be unsorted.
%
% Format: 
% sequence_indicies = fcn_geometry_findPointsInSequence(input_distances, base_point_index, station_tolerance, (fig_num))
%
% INPUTS:
%      input_distances: a Nx1 vector of distances where N is the number of
%      distances.
%
%      base_point_index: the index of the point that anchors the sequence
%
%      station_tolerance: the maximum distance allowed between any points
%      in the sequence either before or after the base point.
%
%      (OPTIONAL INPUTS)
% 
%      fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed.
%
% OUTPUTS:
%
%      sequence_indicies: the indicies that meet the sequence agreement
%      criteria around the base point, including the index of the base
%      point
%
% DEPENDENCIES:
%      (none)
%
% EXAMPLES:
%      
% See the script: script_test_fcn_geometry_findPointsInSequence
% for a full test suite.

% This function was written on 2023_12_29 by S. Brennan
% Questions or comments? sbrennan@psu.edu 

% Revision history:
% 2023_12_29 - S. Brennan
% -- wrote the code

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==4 && isequal(varargin{1},-1))
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
    debug_fig_num = 23434;
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
        narginchk(3,4);

        % Check the input_distances input to be length greater than or
        % equal to 1
        fcn_DebugTools_checkInputsToFunctions(...
            input_distances, '1column_of_numbers', [1 2]);

        % Check the base_point_index input is a positive single integer
        fcn_DebugTools_checkInputsToFunctions(base_point_index, 'positive_1column_of_integers',[1 1]);

        % Check the station_tolerance input is a positive single number
        fcn_DebugTools_checkInputsToFunctions(station_tolerance, 'positive_1column_of_numbers',1);

    end
end

% Does user want to specify fig_num?
fig_num = []; % Default is to have no figure
flag_do_plots = 0;
if 4<= nargin && 0==flag_max_speed
    temp = varargin{end};
    if ~isempty(temp)
        fig_num = temp;
        flag_do_plots = 1;
    end
end


%% Solve for the line fit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _       
%  |  \/  |     (_)      
%  | \  / | __ _ _ _ __  
%  | |\/| |/ _` | | '_ \ 
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Sort the distances in direction of point-to-point projection
[sorted_distances,sort_order_of_indicies] = sort(input_distances);

% sorted_indicies = indicies_in_lateral_agreement(sort_order_of_indicies);
sorted_indicies = sort_order_of_indicies;


% Find the starting point in the sort associated with the base point
% index
starting_point_sort = find(base_point_index == sorted_indicies, 1);
distances_after_base_point = sorted_distances(starting_point_sort:end);
distances_before_base_point = sorted_distances(1:starting_point_sort);
diff_tangent_distances_after_base_point = diff(distances_after_base_point);
diff_tangent_distances_before_base_point = diff(distances_before_base_point);

% Which indicies after the base point are in agreement?
stop_point_after_base_point = find(diff_tangent_distances_after_base_point>station_tolerance,1);
if ~isempty(stop_point_after_base_point)
    sorted_indicies_after_base_point = (starting_point_sort:(starting_point_sort+stop_point_after_base_point-1));
else
    sorted_indicies_after_base_point = (starting_point_sort:length(sorted_indicies));
end

% Which indicies before the base point are in agreement?
stop_point_before_base_point = find(diff_tangent_distances_before_base_point>station_tolerance,1,'last');
if ~isempty(stop_point_before_base_point)
    sorted_indicies_before_base_point = ((stop_point_before_base_point+1):(starting_point_sort-1));
else
    sorted_indicies_before_base_point = (1:(starting_point_sort-1));
end

% Store before and after indicies together to produce indicies in
% station agreement
sequence_sorted_indicies = [sorted_indicies_before_base_point sorted_indicies_after_base_point];
sequence_indicies = sort_order_of_indicies(sequence_sorted_indicies);

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

    hold on;
    grid on;
    title('Distances and sequence match');
    xlabel('Distance [meters]');
    ylabel('Distance [meters]')

    % Plot the input distances 
    plot(input_distances(:,1),input_distances(:,1),'k.','MarkerSize',20);
    plot(input_distances(base_point_index,1),input_distances(base_point_index,1),'b.','MarkerSize',30);

    % Plot the sequence_indicies distances
    plot(input_distances(sequence_indicies,1),input_distances(sequence_indicies,1),'r.-','MarkerSize',15,'LineWidth',3);

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



