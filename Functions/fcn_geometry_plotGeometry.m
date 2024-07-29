function XY_data = fcn_geometry_plotGeometry(plot_type_string, parameters, varargin)
%% fcn_geometry_plotGeometry
% Plots an individual geometry defined by a string name and parameter set.
%
% Format:
%      XY_data = fcn_geometry_plotGeometry(plot_type_string, parameters, (segment_length), (format), (fig_num))
%
% INPUTS:
%
%      plot_type_string: a string indicating the geometry type to plot,
%      such as 'line', 'segment','spiral, or 'arc'. If plot string is
%      empty, or 'none', then nothing is plotted.
%
%      parameters: the parameter set describing the geometry. See
%      fcn_geometry_fillEmptyDomainStructure for details, as the parameter
%      set is different for each geometry type.
%
%      (OPTIONAL INPUTS)
%
%      segment_length: the smallest step to use for plotting, representing
%      the length (approximately) between points. Default is 0.1 meters.
%
%      format: A format string, e.g. 'b-', that dictates the plot style or
%      a color vector, e.g. [1 0 0.23], that dictates the line color. The
%      format string can also be a complex string - see the test script for
%      examples.
%
%      fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed. If left empty, just plots to the current
%      figure.
%
% OUTPUTS:
%
%      XY_data: the data produced during plotting calculations. Note: this
%      data is returned even if fig_num is empty or set to -1.
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
% 2024_04_15 - S. Brennan
% -- added XY_data output
% -- added segment_length input
% 2024_04_20 - S. Brennan
% -- added spiral type
% 2024_04_30 - S. Brennan
% -- checks added to avoid plotting nan values for parameters
% 2024_05_03 - S. Brennan
% -- added line and circle types
% 2024_05_05 - S. Brennan
% -- changed start/end plotting to use dots/circles
% 2024_05_06 - Aneesh Batchu
% -- Added line segment as one of the names for the segment case
% 2024_05_09 - S. Brennan
% -- added 'none' and empty as allowable geometry types
% 2024_05_10 - S. Brennan
% -- added format string input option
% 2024_05_16
% -- Fixed bug that happens when XY data is empty
% 2024_06_16 - Sean Brennan
% -- changed spiral parameter format to new style:
%            'spiral' - 
%
%               [
%                x0,  % The initial x value
%                y0,  % The initial y value
%                h0,  % The initial heading
%                s_Length,  % the s-coordinate length allowed
%                K0,  % The initial curvature
%                Kf   % The final curvature
%              ] 
% 2024_06_19 - Sean Brennan
% -- changed line parameter format to new standard:
%             [
%              base_point_x, 
%              base_point_y, 
%              heading,
%             ]
% 2024_06_19 - Sean Brennan
% -- changed segment parameter format to new standard:
%             [
%              base_point_x, 
%              base_point_y, 
%              heading,
%              s_Length,
%             ]

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==5 && isequal(varargin{end},-1))
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
    if flag_check_inputs == 1
        % Are there the right number of inputs?
        narginchk(2,5);

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

% Does user want to specify the segment_length?
segment_length = 0.1;
if 3 <= nargin
    temp = varargin{1};
    if ~isempty(temp)
        segment_length = temp;
    end
end

% Does user want to change the plot style?
% Set plotting defaults
plot_str = '';
plot_type = 1;  % Plot type refers to 1: a string is given or 2: a color is given - default is 1

% Check to see if user passed in a string or color style?
if 4 <= nargin
    input = varargin{2};
    if ~isempty(input)
        plot_str = input;
        if isnumeric(plot_str)  % Numbers are a color style
            plot_type = 2;
        end
    end
end

% Does user want to specify fig_num?
flag_do_plots = 0;
if 0==flag_max_speed
    flag_do_plots = 1;
    if 5<= nargin
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
% Calculate the XY_data depending on type
if isempty(plot_type_string)
    XY_data = [nan nan; nan nan];
else
    switch lower(plot_type_string)
        case {'none'}
            XY_data = [nan nan; nan nan];
        case {'line'}
            if ~isempty(parameters) && ~any(isnan(parameters))

                base_point_xy        = parameters(1,1:2);
                line_vector          = [cos(parameters(1,3))   sin(parameters(1,3))  ];

                station_distance_min = -10;
                station_distance_max = 10;

                stations = (station_distance_min:segment_length:station_distance_max)';

                if stations(end)~=station_distance_max
                    stations = [stations; station_distance_max];
                end

                XY_data = stations*line_vector + ones(length(stations),1)* base_point_xy;
            else
                XY_data = [nan nan; nan nan];
            end

        case {'segment','vector regression segment fit', 'line segment'}
            if ~isempty(parameters) && ~any(isnan(parameters))
                base_point_xy        = parameters(1,1:2);
                line_vector          = [cos(parameters(1,3))   sin(parameters(1,3))  ];
                station_length       = parameters(1,4);

                if station_length>=0
                    stations = (0:segment_length:station_length)';
                    station_end = station_length;
                else
                    stations = (station_length:segment_length:0)';
                    station_end = 0;
                end
                if stations(end)~=station_end
                    stations = [stations; station_end];
                end

                XY_data = stations*line_vector + ones(length(stations),1)* base_point_xy;
            else
                XY_data = [nan nan; nan nan];
            end

        case {'circle'}
            if ~isempty(parameters) && ~any(isnan(parameters))
                circleCenter         = parameters(1,1:2);
                circleRadius         = parameters(1,3);
                arcAngles            = [0 2*pi];

                start_angle_in_radians = arcAngles(1);
                end_angle_in_radians   = arcAngles(2);
                degree_step = (segment_length/circleRadius)*180/pi;
                flag_arc_is_counterclockwise = 1;

                XY_data = fcn_geometry_plotArc(circleCenter, circleRadius, start_angle_in_radians, end_angle_in_radians, (flag_arc_is_counterclockwise), (degree_step),[],[]);
            else
                XY_data = [nan nan; nan nan];
            end

        case {'arc','regression arc'}
            if ~isempty(parameters) && ~any(isnan(parameters))
                circleCenter         = parameters(1,1:2);
                circleRadius         = parameters(1,3);
                arcAngles            = parameters(1,4:5);
                flag_arc_is_counterclockwise = parameters(1,7);

                start_angle_in_radians = arcAngles(1);
                end_angle_in_radians   = arcAngles(2);
                degree_step = (segment_length/circleRadius)*180/pi;

                XY_data = fcn_geometry_plotArc(circleCenter, circleRadius, start_angle_in_radians, end_angle_in_radians, (flag_arc_is_counterclockwise), (degree_step),[],[]);
            else
                XY_data = [nan nan; nan nan];
            end
        case {'spiral'}
            % 'spiral' -
            % 
            % [
            %     x0,  % The initial x value
            %     y0,  % The initial y value
            %     h0,  % The initial heading
            %     s_Length,  % the s-coordinate length allowed
            %     K0,  % The initial curvature
            %     Kf   % The final curvature
            %     ]
            if ~isempty(parameters) && ~any(isnan(parameters))
                x0                = parameters(1,1);
                y0                = parameters(1,2);
                h0                = parameters(1,3);
                spiralLength      = parameters(1,4);
                K0                = parameters(1,5);
                Kf                = parameters(1,6);

                s_range = linspace(0,spiralLength,ceil(spiralLength/segment_length));
                [x_spiral,y_spiral] = fcn_geometry_extractXYfromSTSpiral(s_range,spiralLength,h0,x0,y0,K0,Kf);

                XY_data = [x_spiral,y_spiral];
            else
                XY_data = [nan nan; nan nan];
            end

        otherwise
            warning('on','backtrace');
            warning('An error will now be thrown because a geometry string was not recognized.');
            error('Unknown plotting type: %s', plot_type_string);
    end
end % Ends check if string is empty


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

    % Check to see if need to ammend color specification to a complex plot_str
    if length(plot_str)>3 && ~contains(lower(plot_str),'color')
        plot_str = cat(2,plot_str,',''Color'',color_vector');
    end

    % Plot lines as quiver arrows
    if strcmp(plot_type_string,'line')
        direction_vector = XY_data(end,:)-XY_data(1,:);
        if plot_type==1
            if length(plot_str)>3
                eval_string = sprintf('quiver(XY_data(1,1),XY_data(1,2),    direction_vector(1,1),direction_vector(1,2), 0, %s)',plot_str);
                eval(eval_string);
                eval_string = sprintf('quiver(XY_data(end,1),XY_data(end,2),-direction_vector(1,1),-direction_vector(1,2), 0, %s)',plot_str);
                eval(eval_string);
            elseif ~isempty(plot_str)
                quiver(XY_data(1,1),XY_data(1,2),    direction_vector(1,1),direction_vector(1,2),0,plot_str);
                quiver(XY_data(end,1),XY_data(end,2),-direction_vector(1,1),-direction_vector(1,2),0,plot_str);
            else
                % Plot string is empty - use defaults
                quiver(XY_data(1,1),XY_data(1,2),    direction_vector(1,1),direction_vector(1,2),0,'-','LineWidth',3,'Color',color_vector);
                quiver(XY_data(end,1),XY_data(end,2),-direction_vector(1,1),-direction_vector(1,2),0,'-','LineWidth',3,'Color',color_vector);
            end
        elseif plot_type==2
            quiver(XY_data(1,1),XY_data(1,2),    direction_vector(1,1),direction_vector(1,2),0,'-','LineWidth',3,'Color',plot_str);
            quiver(XY_data(end,1),XY_data(end,2),-direction_vector(1,1),-direction_vector(1,2),0,'-','LineWidth',3,'Color',plot_str);
        end

    else
        % Plot everything else normally as XY data
        if plot_type==1
            if length(plot_str)>3                
                eval_string = sprintf('plot(XY_data(:,1),XY_data(:,2),%s)',plot_str);
                eval(eval_string);
            elseif ~isempty(plot_str)
                plot(XY_data(:,1),XY_data(:,2),plot_str);
            else
                % Plot string is empty - use defaults
                plot(XY_data(:,1),XY_data(:,2),'-','LineWidth',3,'Color',color_vector);
            end
        elseif plot_type==2
            plot(XY_data(:,1),XY_data(:,2),'-','LineWidth',3,'Color',plot_str);
        end

    end

    % Plot green/red headers and tailers?
    if ~any(isnan(XY_data),'all') && ~isempty(XY_data)
        if 1==1
            %%%%%%
            % Plot start and end locations as green/red points
            plot(XY_data(1,1),XY_data(1,2),     '.','Color',[0 1 0],'Linewidth',5,'MarkerSize',20);
            plot(XY_data(end,1),XY_data(end,2), 'o','Color',[1 0 0],'MarkerSize',10);

        else
            %%%%%%
            % Plot as start and end locations as green/red line segments

            % Plot green headers - calculated from vector direction
            maximum_arrow_length = 2*segment_length;
            minimum_arrow_length = 1*segment_length;

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
        end
    end


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


