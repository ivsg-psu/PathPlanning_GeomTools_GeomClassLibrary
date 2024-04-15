function fcn_geometry_plotGeometry(plot_type_string, parameters, varargin)
%% fcn_geometry_plotGeometry
% Plots an individual geometry defined by a string name and parameter set.
%
% Format:
% fcn_geometry_plotGeometry(plot_type_string, parameters, (fig_num))
%
% INPUTS:
%
%      plot_type_string: a string indicating the geometry type to plot,
%      such as 'line' or 'arc'.
%
%      parameters: the parameter set describing the geometry. See
%      fcn_geometry_fillEmptyDomainStructure for details, as the parameter
%      set is different for each geometry type.
%
%      (OPTIONAL INPUTS)
%
%      fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed. If left empty, just plots to the current
%      figure.
%
% OUTPUTS:
%
%      (none)
%
% DEPENDENCIES:
%      
%      (none)
%
% EXAMPLES:
%
% See the script: script_test_fcn_geometry_plotGeometry
% for a full test suite.
%
% This function was written on 2024_04_14 by S. Brennan
% Questions or comments? sbrennan@psu.edu

% Revision history:
% 2024_04_14 - S. Brennan
% -- wrote the code

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==3 && isequal(varargin{end},-1))
    flag_do_debug = 0; % % Flag to plot the results for debugging
    flag_check_inputs = 0; % Flag to perform input checking
    flag_max_speed = 1;
else
    % Check to see if we are externally setting debug mode to be "on"
    flag_do_debug = 0; % % Flag to plot the results for debugging
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
    debug_fig_num = 34838; 
else
    debug_fig_num = []; 
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
        narginchk(2,3);

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
if 0==flag_max_speed
    flag_do_plots = 1;
    if 3<= nargin
        temp = varargin{end};
        if ~isempty(temp)
            fig_num = temp;
        end
    else
        fig_num = gcf;
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
% Nothing in main, all is in plotting


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
    xlabel('X [m]');
    ylabel('Y [m]');

    % Get the color vector using the name
    color_vector = fcn_geometry_fillColorFromNumberOrName(2,lower(plot_type_string));

    % Change the plot style depending on type
    switch lower(plot_type_string)
        case {'line','segment','vector regression segment fit'}
            line_vector          = parameters(1,1:2);
            base_point_xy        = parameters(1,3:4);
            station_distance_min = parameters(1,5);
            station_distance_max = parameters(1,6);

            stations = [station_distance_min; station_distance_max];
            XY_data = stations*line_vector + ones(length(stations),1)* base_point_xy;
            plot(XY_data(:,1),XY_data(:,2),'-','LineWidth',3,'Color',color_vector);

        case {'arc','regression arc'}
            circleCenter         = parameters(1,1:2);
            circleRadius         = parameters(1,3);
            arcAngles            = parameters(1,4:5);
            flag_arc_is_counterclockwise = parameters(1,7);

            start_angle_in_radians = arcAngles(1);
            end_angle_in_radians   = arcAngles(2);
            degree_step = [];

            XY_data = fcn_geometry_plotArc(circleCenter, circleRadius, start_angle_in_radians, end_angle_in_radians, (flag_arc_is_counterclockwise), (degree_step),[],[]);


            plot(XY_data(:,1),XY_data(:,2),'-','LineWidth',3,'Color',color_vector);

        otherwise
            warning('on','backtrace');
            warning('An error will now be thrown because a geometry string was not recognized.');
            error('Unknown plotting type: %s', plot_type_string);
    end

    % Plot green headers - calculated from vector direction

    maximum_arrow_length = 10;
    minimum_arrow_length = 0.2;

    vector_direction_start = XY_data(2,1:2) - XY_data(1,1:2);
    start_length = sum(vector_direction_start.^2,2).^0.5;
    unit_vector_direction_start = fcn_geometry_calcUnitVector(vector_direction_start);
    arrow_length = max(min(maximum_arrow_length,start_length*0.2),minimum_arrow_length);
    offset_start = XY_data(1,1:2) + arrow_length*unit_vector_direction_start;

    start_line = [XY_data(1,1:2) 0; offset_start, 0];
    plot(start_line(:,1),start_line(:,2), '-','Color',[0 1 0],'Linewidth',5);

    % Plot red tailers - calculated from vector direction
    vector_direction_end = (XY_data(end,1:2) - XY_data(end-1,1:2));
    end_length = sum(vector_direction_end.^2,2).^0.5;
    unit_vector_direction_end = fcn_geometry_calcUnitVector(vector_direction_end);
    arrow_length = max(min(maximum_arrow_length,end_length*0.2),minimum_arrow_length);
    offset_end = XY_data(end,1:2) - arrow_length*unit_vector_direction_end;
    end_line = [offset_end, 0; XY_data(end,1:2) 0];
    plot(end_line(:,1),end_line(:,2), '-','Color',[1 0 0],'Linewidth',5);



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


