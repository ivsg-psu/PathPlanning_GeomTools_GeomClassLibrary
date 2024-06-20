function [XY_parameters] = ...
fcn_geometry_orientGeometrySt2XY(st_parameters_type_strings, st_parameters, St_transform_XYtoSt, flag_primary_parameter_is_flipped, varargin)
%% fcn_geometry_orientGeometrySt2XY
% Rotates geometries from XY into ST coordinates so that the ending
% geometry of a primary parameter set is tangent to the X-axis and ends at
% the origin. 
%
% This is generally used as a transform prior to geometric manipulations to
% ease complexity and permutations needed to test geometric manipulation
% codes.
%
% Format:
% [xy_parameters] = ...
% fcn_geometry_orientGeometrySt2XY(st_parameters_type_strings, st_parameters, St_transform, flag_primary_parameter_is_flipped, (fig_num))
%
% INPUTS:
%
%      st_parameters_type_strings: a cell array of N strings containing the
%      parameter type of the parameters to be returned to xy orientation
%      from st. The string can be one of the basic geometries defined in
%      fcn_geometry_fillEmptyDomainStructure:
% 
%           'line': this is a line.
%
%           'segment': this is a line segment.
%
%           'circle': this is a circle 
%
%           'arc': this is an arc 
%
%           'spiral': this is a spiral 
%
%      st_parameters: a cell array of the N parameter vectors describing
%      the st-geometry for each of the geometries gifen by the
%      st_parameter_type_strings.
%
%      St_transform_XYtoSt: the 4x4 transformation matrix compatible with
%      the SE library and 'transform' function usage that describes the
%      transformation applied to each geometry to produce St from XY. Note
%      that this is produced in the XY to St transform.
%
%      flag_primary_parameter_is_flipped: for an arc geometry, if the arc
%      is clock-wise, the arc will be flipped about the x-axis such that it
%      is always counter-clockwise, and all secondary geometries will also
%      be flipped. If this "flipped" situation was used to create the St
%      geometry, flag_primary_parameter is set to 1. Otherwise, it is 0.
%
%      (OPTIONAL INPUTS)
%
%      fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed.
%
% OUTPUTS:
%
%      xt_parameters: a cell array of each of the parameter sets
%      describing the St geometries in XY coordinates
%
% DEPENDENCIES:
%
%      (none)
%
% EXAMPLES:
%
% See the script: script_test_fcn_geometry_orientGeometrySt2XY
% for a full test suite.
%
% This function was written on 2024_05_04 by S. Brennan
% Questions or comments? sbrennan@psu.edu

% Revision history:
% 2024_05_04 - S. Brennan
% -- wrote the code
% 2024_05_09 - S. Brennan
% -- fixed bug in segment calculation wherein unit vector gives NaN if
% start and end points are same
% -- fixed bug where NaN inputs cause it to crash
% 2024_06_16 - Sean Brennan
% -- changed parameter format to new style:
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
% 2024_06_16 - Sean Brennan
% -- changed parameter format to new style:
%            'spiral' - 
%               [
%                x0,  % The initial x value
%                y0,  % The initial y value
%                h0,  % The initial heading
%                s_Length,  % the s-coordinate length allowed
%                K0,  % The initial curvature
%                Kf   % The final curvature
%              ] 
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
    flag_do_debug = 0; % % % Flag to plot the results for debugging
    flag_check_inputs = 0; % Flag to perform input checking
    flag_max_speed = 1;
else
    % Check to see if we are externally setting debug mode to be "on"
    flag_do_debug = 0; % % % Flag to plot the results for debugging
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
        narginchk(4,5);

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
if 5<= nargin && 0==flag_max_speed
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
% 
% flag_primary_parameter_is_flipped = 0;
% switch st_parameters_type_strings
%     case 'segment'
%         line_unit_tangent_vector   = st_parameters(1,1:2);
%         line_base_point_xy         = st_parameters(1,3:4);
%         % line_s_start               = primary_parameters(1,5);
%         line_s_end                 = st_parameters(1,6);
%         % line_start_xy              = line_base_point_xy + line_unit_tangent_vector*line_s_start;
%         line_end_xy                = line_base_point_xy + line_unit_tangent_vector*line_s_end;
% 
%         % Use the SE2 toolbox to transform 
%         % Start with translation of the line's center to the origin 
%         translation_to_origin        = -line_end_xy;          % Push end of the segment to be at origin
%         transformMatrix_translation_into_origin = se2(0,'theta',translation_to_origin);
% 
%         % Rotate everything to align the line to x-axis
%         line_angle = atan2(line_unit_tangent_vector(2),line_unit_tangent_vector(1));
%         rotation_angle = -1*(line_angle);
%         transformMatrix_rotation_into_St         = se2(rotation_angle,'theta',[0 0]);
% 
%         % Combine all the transformations
%         St_transform  = transformMatrix_rotation_into_St*transformMatrix_translation_into_origin;
% 
%     case 'arc'
% 
%         arc_center_xy                = st_parameters(1,1:2);
%         arc_radius                   = st_parameters(1,3);       
%         arc_end_angle_in_radians     = st_parameters(1,5);
% 
% 
%         arc_is_counter_clockwise     = st_parameters(1,7);
% 
%         if arc_is_counter_clockwise==1
%             counter_clockwise_multiplier = 1;
%         else
%             counter_clockwise_multiplier = -1;
%         end
% 
%         % Use the SE2 toolbox to transform 
%         % Start with translation of the arc's center to the origin 
%         translation_to_center1        = -arc_center_xy;          % Push circle 1 to be at origin
%         transformMatrix_translation_into_center1 = se2(0,'theta',translation_to_center1);
% 
%         % Rotate everything to push the end of the arc1 angle to -90
%         % degrees
%         rotation_angle = -1*(arc_end_angle_in_radians + counter_clockwise_multiplier*pi/2);
%         transformMatrix_rotation_into_St         = se2(rotation_angle,'theta',[0 0]);
% 
%         % Finally, translate everything upwards (counter-clockwise) or downwards (clockwise) to make arc1's end at the
%         % origin
%         translation_to_offset_center1  = [0 counter_clockwise_multiplier*arc_radius];   
%         transformMatrix_offset_center1 = se2(0,'theta',translation_to_offset_center1);
% 
%         % Is the arc flipped about y-axis?
%         transformMatrix_rotation_about_xaxis = se2(eye(3));
%         if arc_is_counter_clockwise~=1
%             transformMatrix_rotation_about_xaxis = se2(diag([1 -1 1],0));
%             flag_primary_parameter_is_flipped = 1;
%         end
% 
%         % Combine all the transformations
%         St_transform  = transformMatrix_rotation_about_xaxis*transformMatrix_offset_center1*transformMatrix_rotation_into_St*transformMatrix_translation_into_center1;        
% 
%     otherwise
%         warning('on','backtrace');
%         warning('An error will be thrown at this point due to incorrect primary_parameters_type_string.');
% 
%         error('ST alignments are not yet supported for geometries of type: %s',st_parameters_type_strings);
% end


%% Convert all the fits based on the rotations and transformations

% Check to see if the input array is a cell type
if ~iscell(st_parameters_type_strings) 
    if (isstring(st_parameters_type_strings) || ischar(st_parameters_type_strings))
        temp = st_parameters_type_strings;
        clear st_parameters_type_strings;
        st_parameters_type_strings{1} = temp;
    else
        error('Unrecongized type for st_parameters_type_strings. Quitting.');
    end
end
if ~iscell(st_parameters) 
    if isnumeric(st_parameters)
        temp = st_parameters;
        clear st_parameters;
        st_parameters{1} = temp;
    else
        error('Unrecongized type for st_parameters. Quitting.');
    end
end





N_conversions = length(st_parameters_type_strings);

% Initialize the outputs and transforms
XY_parameters      = st_parameters;
transform_St2XY = inv(St_transform_XYtoSt);

% Loop through all the parameters, converting each depending on their
% respective type
for ith_parameter_set = 1:N_conversions
    
    % Get the current parameter set
    current_st_parameters = st_parameters{ith_parameter_set};

    % Check if parameters contain NaN values
    if any(isnan(current_st_parameters))
        XY_current_parameters = nan(size(current_st_parameters));
    else
        % Parameters do not contain NaN, so convert them
        current_parameters_type = st_parameters_type_strings{ith_parameter_set};

        switch current_parameters_type
            case 'line'

                % Get the line details from parameters
                line_base_point_xy         = current_st_parameters(1,1:2);
                line_unit_tangent_vector   = [cos(current_st_parameters(1,3)) sin(current_st_parameters(1,3))];

                %%%%%
                % Fix line
                line_unit_tangent_vector_end    = line_base_point_xy+line_unit_tangent_vector;
                line_unit_tangent_vector_end_St = transform(transform_St2XY,line_unit_tangent_vector_end);
                line_base_point_St              = transform(transform_St2XY,line_base_point_xy);

                st_line_parameters(1,1:2)       = line_base_point_St;
                line_unit_vector                = fcn_geometry_calcUnitVector(line_unit_tangent_vector_end_St - line_base_point_St);
                st_line_parameters(1,3)         = atan2(line_unit_vector(2),line_unit_vector(1));

                XY_current_parameters = st_line_parameters;

            case 'segment'
                % Get the segment details from parameters
                segment_base_point_xy       = current_st_parameters(1,1:2);
                segment_unit_tangent_vector = [cos(current_st_parameters(1,3)) sin(current_st_parameters(1,3))];
                segment_length              = current_st_parameters(1,4);

                %%%%%
                % Fix segment
                segment_start_point_xy        = segment_base_point_xy; 
                segment_start_point_St        = transform(transform_St2XY,segment_start_point_xy);
                segment_unitVector_start_point_xy        = [0 0];
                segment_unitVector_start_point_St        = transform(transform_St2XY,segment_unitVector_start_point_xy);
                segment_unitVector_end_point_xy          = segment_unit_tangent_vector;
                segment_unitVector_end_point_St          = transform(transform_St2XY,segment_unitVector_end_point_xy);

                st_segment_parameters(1,1:2)  = segment_start_point_St;
                segment_unit_vector = fcn_geometry_calcUnitVector(segment_unitVector_end_point_St - segment_unitVector_start_point_St);
                st_segment_parameters(1,3)    = atan2(segment_unit_vector(2),segment_unit_vector(1));
                st_segment_parameters(1,4)    = segment_length;

                XY_current_parameters = st_segment_parameters;


            case 'circle'
                % A test circle
                % [circleCenter_x.
                %     circleCenter_y,
                %     radius]

                % Get the circle details from parameters
                circle_center_xy                = current_st_parameters(1,1:2);
                circle_radius                   = current_st_parameters(1,3);

                %%%%%
                % Fix circle
                circle_center_St                 = transform(transform_St2XY,circle_center_xy);


                circle_parameters(1,1:2)      = circle_center_St;
                circle_parameters(1,3)        = circle_radius;

                XY_current_parameters = circle_parameters;


            case 'arc'
                % Get the arc details from parameters
                arc_center_st                = current_st_parameters(1,1:2);
                arc_radius                   = current_st_parameters(1,3);
                arc_start_angle_in_radians   = current_st_parameters(1,4);
                arc_end_angle_in_radians     = current_st_parameters(1,5);
                arc_is_circle                = current_st_parameters(1,6);
                arc_is_counter_clockwise     = current_st_parameters(1,7);

                %%%%%
                % Fix arc
                arc_center_XY                 = transform(transform_St2XY,arc_center_st);
                arc_start_st                  = arc_center_st + arc_radius*[cos(arc_start_angle_in_radians) sin(arc_start_angle_in_radians)];
                arc_end_st                    = arc_center_st + arc_radius*[cos(arc_end_angle_in_radians)   sin(arc_end_angle_in_radians)];
                arc_start_XY                  = transform(transform_St2XY,arc_start_st);
                arc_end_XY                    = transform(transform_St2XY,arc_end_st);
                arc_start_vector_XY           = arc_start_XY - arc_center_XY;
                arc_end_vector_XY             = arc_end_XY   - arc_center_XY;


                st_arc_parameters(1,1:2)      = arc_center_XY;
                st_arc_parameters(1,3)        = arc_radius;
                st_arc_parameters(1,4)        = atan2(arc_start_vector_XY(2),arc_start_vector_XY(1));
                st_arc_parameters(1,5)        = atan2(arc_end_vector_XY(2),  arc_end_vector_XY(1));
                st_arc_parameters(1,6)        = arc_is_circle;
                if  flag_primary_parameter_is_flipped ~= 1
                    st_arc_parameters(1,7)        = arc_is_counter_clockwise;
                else
                    if arc_is_counter_clockwise
                        st_arc_parameters(1,7)        = 0;
                    else
                        st_arc_parameters(1,7)        = 1;
                    end
                end

                XY_current_parameters = st_arc_parameters;

            case 'spiral'
                % A test spiral
                %  [
                %   x0,  % The initial x value
                %   y0,  % The initial y value
                %   h0,  % The initial heading
                %   spiralLength,  % the s-coordinate length allowed
                %   K0,  % The initial curvature
                %   Kf   % The final curvature
                % ]

                % Get the spiral details from parameters
                spiral_center_ST             = current_st_parameters(1,1:2);
                spiral_heading               = current_st_parameters(1,3);
                spiral_length                = current_st_parameters(1,4);
                spiral_K0                    = current_st_parameters(1,5);
                spiral_Kf                    = current_st_parameters(1,6);

                %%%%%
                % Fix spiral
                spiral_heading_vector_root_ST = [0 0];
                spiral_heading_vector_root_XY = transform(transform_St2XY,spiral_heading_vector_root_ST);
                spiral_heading_vector_head_ST = [cos(spiral_heading) sin(spiral_heading)];
                spiral_heading_vector_head_XY = transform(transform_St2XY,spiral_heading_vector_head_ST);
                spiral_heading_vector_XY      = spiral_heading_vector_head_XY - spiral_heading_vector_root_XY;
                spiral_heading_XY             = atan2(spiral_heading_vector_XY(2),spiral_heading_vector_XY(1));

                spiral_center_XY              = transform(transform_St2XY,spiral_center_ST);

                if  flag_primary_parameter_is_flipped ~= 1
                    spiral_K0_XY    = spiral_K0;
                    spiral_Kf_XY    = spiral_Kf;
                else
                    spiral_K0_XY    = -spiral_K0;
                    spiral_Kf_XY    = -spiral_Kf;
                end


                spiral_parameters(1,1:2)      = spiral_center_XY;
                spiral_parameters(1,3)        = spiral_heading_XY;
                spiral_parameters(1,4)        = spiral_length;
                spiral_parameters(1,5)        = spiral_K0_XY;
                spiral_parameters(1,6)        = spiral_Kf_XY;

                XY_current_parameters = spiral_parameters;


            otherwise
                warning('on','backtrace');
                warning('An error will be thrown at this point due to unknown parameter type.');
                error('Alignments are not yet supported for curves from fit type: %s',current_parameters_type);
        end
    end % Ends check to see if parameters do not contain NaN
    XY_parameters{ith_parameter_set} = XY_current_parameters;
end % Ends looping through the parameter sets.



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

    subplot(1,2,1);
    hold on;
    grid on;
    axis equal;
    sgtitle('Converting ST to XY');
    title('ST');




    % Loop through all the parameters
    for ith_parameter_set = 1:N_conversions

        % Get the current parameter set
        current_parameters_type = st_parameters_type_strings{ith_parameter_set};
        current_st_parameters = st_parameters{ith_parameter_set};
       
        % Plot the inputs
        fcn_geometry_plotGeometry(current_parameters_type,current_st_parameters);
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
    xlabel('S [meters]');
    ylabel('t [meters]')

    temp_axis = axis;

    subplot(1,2,2);
    hold on;
    grid on;
    axis equal;
    title('XY');

    % Plot the outputs
    % Loop through all the parameters
    for ith_parameter_set = 1:N_conversions

        % Get the current parameter set
        current_parameters_type = st_parameters_type_strings{ith_parameter_set};
        current_st_parameters = XY_parameters{ith_parameter_set};
       
        % Plot the inputs
        fcn_geometry_plotGeometry(current_parameters_type,current_st_parameters);
    end

    % Make axis slightly larger?
    if flag_rescale_axis
        temp = axis;
        %     temp = [min(points(:,1)) max(points(:,1)) min(points(:,2)) max(points(:,2))];
        axis_range_x = temp(2)-temp(1);
        axis_range_y = temp(4)-temp(3);
        percent_larger = 0.3;
        temp_axis2 = [temp(1)-percent_larger*axis_range_x, temp(2)+percent_larger*axis_range_x,  temp(3)-percent_larger*axis_range_y, temp(4)+percent_larger*axis_range_y];
        good_axis = [min([temp_axis(1) temp_axis2(1)])  max([temp_axis(2) temp_axis2(2)])  min([temp_axis(3) temp_axis2(3)])  max([temp_axis(4) temp_axis2(4)]) ];
        axis(good_axis);
    end
    xlabel('X [meters]');
    ylabel('Y [meters]')

    % Force the subplots to have matching axes
    subplot(1,2,2);
    good_axis = axis;

    subplot(1,2,1);
    axis(good_axis)

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

