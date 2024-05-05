function [st_primary_parameters, st_secondary_parameters, St_transform, rotation_angle, flag_primary_parameter_is_flipped] = ...
fcn_orientGeometryXY2St(primary_parameters_type_string, primary_parameters, varargin)
%% fcn_orientGeometryXY2St
% Rotates geometries from XY into ST coordinates so that the ending
% geometry of a primary parameter set is tangent to the X-axis and ends at
% the origin. 
%
% This is generally used as a transform prior to geometric manipulations to
% ease complexity and permutations needed to test geometric manipulation
% codes.
%
% Format:
% [st_primary_parameters, st_secondary_parameters, St_transform, rotation_angle, flag_primary_parameter_is_flipped] = ...
% fcn_orientGeometryXY2St(primary_parameters_type_string, primary_parameters, (secondary_parameters_type_strings), (secondary_parameters), (fig_num))
%
% INPUTS:
%
%      primary_parameters_type_string: a string containing the parameter
%      type of the primary parameters used for orientation. The string can
%      be one of:
% 
%           'segment': this is a line segment that causes a transform to be
%           generated such that the end of the segment will be at the
%           origin, and the segment is on the x-axis.
%
%           'arc': this is an arc that causes a transformation to be
%           generated such that the end of the arc will be at the origin,
%           with the projection of the arc at the end-point being along the
%           x-axis.
%
%      primary_parameters: the parameter set describing the primary
%      geometry. See fcn_geometry_fillEmptyDomainStructure for details for
%      the meaning of each parameter geometry vector.
%
%      (OPTIONAL INPUTS)
%
%      secondary_parameters_type_strings: a cell array of N strings, one
%      string for each type of secondary geometry type. The strings can be
%      one of the following: 'line','segment','circle','arc','spiral'. Each
%      of these geometries will be transformed under the same transform
%      generated from the primary parameter set.
%
%      secondary_parameters: the parameter sets as a cell array of N
%      vectors describing the geometry for each secondary parameter set
%      listed in the secondary_parameters_type_strings. See
%      fcn_geometry_fillEmptyDomainStructure for details.
%
%      fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed.
%
% OUTPUTS:
%
%      st_primary_parameters: the parameter set describing the primary
%      geometry in S-t coordinates.
%
%      st_secondary_parameters: a cell array of each of the parameter sets
%      describing the secondary geometries in the S-t coordinates under the
%      same transformation as the primary_parameters.
%
%      St_transform: the 4x4 transformation matrix compatible with the SE
%      library and 'transform' function usage that describes the
%      transformation applied to each geometry. Note that the inverse
%      transformation is obtained by inv(St_tranform).
%
%      rotation_angle: the angle, in radians, of the rotational component of the transformation.
% 
%      flag_primary_parameter_is_flipped: for an arc geometry, if the arc
%      is clock-wise, the arc will be flipped about the x-axis such that it
%      is always counter-clockwise, and all secondary geometries will also
%      be flipped. If this "flipped" situation is encountered,
%      flag_primary_parameter is set to 1. Otherwise, it is 0.
%
% DEPENDENCIES:
%
%      (none)
%
% EXAMPLES:
%
% See the script: script_test_fcn_orientGeometryXY2St
% for a full test suite.
%
% This function was written on 2024_05_02 by S. Brennan
% Questions or comments? sbrennan@psu.edu

% Revision history:
% 2024_05_02 - S. Brennan
% -- wrote the code

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

% Does user want to specify secondary_parameters_type_strings?
secondary_parameters_type_strings = {};
if (3<=nargin)
    temp = varargin{1};
    if ~isempty(temp)
        secondary_parameters_type_strings = temp;
    end
end


% Does user want to specify secondary_parameters_type_strings?
secondary_parameters = {};
if (4<=nargin)
    temp = varargin{2};
    if ~isempty(temp)
        secondary_parameters = temp;
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

flag_primary_parameter_is_flipped = 0;
switch primary_parameters_type_string
    case 'segment'
        line_unit_tangent_vector   = primary_parameters(1,1:2);
        line_base_point_xy         = primary_parameters(1,3:4);
        % line_s_start               = primary_parameters(1,5);
        line_s_end                 = primary_parameters(1,6);
        % line_start_xy              = line_base_point_xy + line_unit_tangent_vector*line_s_start;
        line_end_xy                = line_base_point_xy + line_unit_tangent_vector*line_s_end;

        % Use the SE2 toolbox to transform 
        % Start with translation of the line's center to the origin 
        translation_to_origin        = -line_end_xy;          % Push end of the segment to be at origin
        transformMatrix_translation_into_origin = se2(0,'theta',translation_to_origin);

        % Rotate everything to align the line to x-axis
        line_angle = atan2(line_unit_tangent_vector(2),line_unit_tangent_vector(1));
        rotation_angle = -1*(line_angle);
        transformMatrix_rotation_into_St         = se2(rotation_angle,'theta',[0 0]);

        % Combine all the transformations
        St_transform  = transformMatrix_rotation_into_St*transformMatrix_translation_into_origin;

    case 'arc'

        arc1_center_xy                = primary_parameters(1,1:2);
        arc1_radius                   = primary_parameters(1,3);       
        arc1_end_angle_in_radians     = primary_parameters(1,5);

        
        arc1_is_counter_clockwise     = primary_parameters(1,7);

        if arc1_is_counter_clockwise==1
            counter_clockwise_multiplier = 1;
        else
            counter_clockwise_multiplier = -1;
        end

        % Use the SE2 toolbox to transform 
        % Start with translation of the arc's center to the origin 
        translation_to_center1        = -arc1_center_xy;          % Push circle 1 to be at origin
        transformMatrix_translation_into_center1 = se2(0,'theta',translation_to_center1);

        % Rotate everything to push the end of the arc1 angle to -90
        % degrees
        rotation_angle = -1*(arc1_end_angle_in_radians + counter_clockwise_multiplier*pi/2);
        transformMatrix_rotation_into_St         = se2(rotation_angle,'theta',[0 0]);

        % Finally, translate everything upwards (counter-clockwise) or downwards (clockwise) to make arc1's end at the
        % origin
        translation_to_offset_center1  = [0 counter_clockwise_multiplier*arc1_radius];   
        transformMatrix_offset_center1 = se2(0,'theta',translation_to_offset_center1);

        % Is the arc flipped about y-axis?
        transformMatrix_rotation_about_xaxis = se2(eye(3));
        if arc1_is_counter_clockwise~=1
            transformMatrix_rotation_about_xaxis = se2(diag([1 -1 1],0));
            flag_primary_parameter_is_flipped = 1;
        end

        % Combine all the transformations
        St_transform  = transformMatrix_rotation_about_xaxis*transformMatrix_offset_center1*transformMatrix_rotation_into_St*transformMatrix_translation_into_center1;        

    otherwise
        warning('on','backtrace');
        warning('An error will be thrown at this point due to incorrect primary_parameters_type_string.');

        error('ST alignments are not yet supported for geometries of type: %s',primary_parameters_type_string);
end


%% Convert all the fits based on the rotations and transformations
N_conversions = 1 + length(secondary_parameters_type_strings);

% Fill in all the types we need to convert, starting with the primary
% parameters in the 1st spot
parameters_to_convert_type_strings{N_conversions} = [];
parameters_to_convert{N_conversions}              = [];
parameters_to_convert_type_strings{1}             = primary_parameters_type_string;
parameters_to_convert{1}                          = primary_parameters;
St_parameters                                     = parameters_to_convert;

% Fill in the rest of the parameters to convert
for ith_parameter_set = 2:N_conversions
    parameters_to_convert_type_strings{ith_parameter_set} = secondary_parameters_type_strings{ith_parameter_set-1};
    parameters_to_convert{ith_parameter_set} = secondary_parameters{ith_parameter_set-1};
end

% Loop through all the parameters, converting each depending on their
% respective type
for ith_parameter_set = 1:N_conversions
    
    % Get the current parameter set
    current_parameters = parameters_to_convert{ith_parameter_set};
    current_parameters_type = parameters_to_convert_type_strings{ith_parameter_set};

    switch current_parameters_type
        case 'line'
            % A test line
            % [unit_projection_vector_x,
            %     unit_projection_vector_y,
            %     base_point_x,
            %     base_point_y,
            %     ]

            % Get the line details from parameters
            line_unit_tangent_vector   = current_parameters(1,1:2);
            line_base_point_xy         = current_parameters(1,3:4);

            %%%%%
            % Fix line
            line_unit_tangent_vector_end    = line_base_point_xy+line_unit_tangent_vector;
            line_unit_tangent_vector_end_St = transform(St_transform,line_unit_tangent_vector_end);
            line_base_point_St              = transform(St_transform,line_base_point_xy);

            st_line_parameters(1,1:2)       = line_unit_tangent_vector_end_St - line_base_point_St;
            st_line_parameters(1,3:4)       = line_base_point_St;

            St_current_parameters = st_line_parameters;

        case 'segment'
            % Get the segment details from parameters
            segment_unit_tangent_vector   = current_parameters(1,1:2);
            segment_base_point_xy         = current_parameters(1,3:4);
            segment_s_start               = current_parameters(1,5);
            segment_s_end                 = current_parameters(1,6);
            segment_length                = segment_s_end - segment_s_start;

            %%%%%
            % Fix segment
            segment_start_point_xy        = segment_base_point_xy + segment_unit_tangent_vector*segment_s_start;
            segment_start_point_St        = transform(St_transform,segment_start_point_xy);
            segment_end_point_xy          = segment_base_point_xy + segment_unit_tangent_vector*segment_s_end;
            segment_end_point_St          = transform(St_transform,segment_end_point_xy);

            st_segment_parameters(1,1:2)  = fcn_geometry_calcUnitVector(segment_end_point_St - segment_start_point_St);
            st_segment_parameters(1,3:4)  = segment_start_point_St;
            st_segment_parameters(1,5)    = 0;
            st_segment_parameters(1,6)    = segment_length;

            St_current_parameters = st_segment_parameters;

        case 'circle'
            % A test circle
            % [circleCenter_x.
            %     circleCenter_y,
            %     radius]

            % Get the circle details from parameters
            circle_center_xy                = current_parameters(1,1:2);
            circle_radius                   = current_parameters(1,3);
            
            %%%%%
            % Fix circle
            circle_center_St                 = transform(St_transform,circle_center_xy);


            circle_parameters(1,1:2)      = circle_center_St;
            circle_parameters(1,3)        = circle_radius;

            St_current_parameters = circle_parameters;

       
        case 'arc'
            % Get the arc details from parameters
            arc1_center_xy                = current_parameters(1,1:2);
            arc1_radius                   = current_parameters(1,3);
            arc1_start_angle_in_radians   = current_parameters(1,4);
            arc1_end_angle_in_radians     = current_parameters(1,5);
            arc1_is_circle                = current_parameters(1,6);
            arc1_is_counter_clockwise     = current_parameters(1,7);

            %%%%%
            % Fix arc
            arc1_center_St                 = transform(St_transform,arc1_center_xy);
            arc1_start_xy                  = arc1_center_xy + arc1_radius*[cos(arc1_start_angle_in_radians) sin(arc1_start_angle_in_radians)];
            arc1_end_xy                    = arc1_center_xy + arc1_radius*[cos(arc1_end_angle_in_radians)   sin(arc1_end_angle_in_radians)];
            arc1_start_St                  = transform(St_transform,arc1_start_xy);
            arc1_end_St                    = transform(St_transform,arc1_end_xy);
            arc1_start_vector_St           = arc1_start_St - arc1_center_St;
            arc1_end_vector_St             = arc1_end_St   - arc1_center_St;


            st_arc1_parameters(1,1:2)      = arc1_center_St;
            st_arc1_parameters(1,3)        = arc1_radius;
            st_arc1_parameters(1,4)        = atan2(arc1_start_vector_St(2),arc1_start_vector_St(1));
            st_arc1_parameters(1,5)        = atan2(arc1_end_vector_St(2),  arc1_end_vector_St(1));
            st_arc1_parameters(1,6)        = arc1_is_circle;
            if  flag_primary_parameter_is_flipped ~= 1
                st_arc1_parameters(1,7)        = arc1_is_counter_clockwise;
            else
                if arc1_is_counter_clockwise
                    st_arc1_parameters(1,7)        = 0;
                else
                    st_arc1_parameters(1,7)        = 1;
                end
            end

            St_current_parameters = st_arc1_parameters;

        case 'spiral'
            % A test spiral
            %  [spiralLength,  % the s-coordinate length allowed
            %   h0,  % The initial heading
            %   x0,  % The initial x value
            %   y0,  % The initial y value
            %   K0,  % The initial curvature
            %   Kf   % The final curvature
            % ]

            % Get the spiral details from parameters
            spiral_length                = current_parameters(1,1);
            spiral_heading               = current_parameters(1,2);
            spiral_center_xy             = current_parameters(1,3:4);
            spiral_K0                    = current_parameters(1,5);
            spiral_Kf                    = current_parameters(1,6);

            %%%%%
            % Fix spiral
            spiral_heading_vector_root_xy = [0 0];
            spiral_heading_vector_root_St = transform(St_transform,spiral_heading_vector_root_xy);
            spiral_heading_vector_head_xy = [cos(spiral_heading) sin(spiral_heading)];
            spiral_heading_vector_head_St = transform(St_transform,spiral_heading_vector_head_xy);
            spiral_heading_vector_St      = spiral_heading_vector_head_St - spiral_heading_vector_root_St;
            spiral_heading_St             = atan2(spiral_heading_vector_St(2),spiral_heading_vector_St(1));

            spiral_center_St              = transform(St_transform,spiral_center_xy);

            if  flag_primary_parameter_is_flipped ~= 1
                spiral_K0_St    = spiral_K0;
                spiral_Kf_St    = spiral_Kf;
            else
                spiral_K0_St    = -spiral_K0;
                spiral_Kf_St    = -spiral_Kf;
            end


            spiral_parameters(1,1)        = spiral_length;
            spiral_parameters(1,2)        = spiral_heading_St;
            spiral_parameters(1,3:4)      = spiral_center_St;
            spiral_parameters(1,5)        = spiral_K0_St;
            spiral_parameters(1,6)        = spiral_Kf_St;

            St_current_parameters = spiral_parameters;


        otherwise
            warning('on','backtrace');
            warning('An error will be thrown at this point due to unknown parameter type.');
            error('Alignments are not yet supported for curves from fit type: %s',current_parameters_type);
    end

    St_parameters{ith_parameter_set} = St_current_parameters;
end % Ends looping through the parameter sets.

%% Save outputs
st_primary_parameters = St_parameters{1};

if N_conversions == 1
    st_secondary_parameters{1} = {};
else
    st_secondary_parameters{N_conversions-1} = [];
end

% Fill in the rest of the parameters to convert
for ith_parameter_set = 1:N_conversions-1
    st_secondary_parameters{ith_parameter_set} = St_parameters{ith_parameter_set+1}; %#ok<AGROW>
end



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
    sgtitle(sprintf('Converting XY to ST for base parameter: %s',primary_parameters_type_string));
    title('Inputs');
    xlabel('X [meters]');
    ylabel('Y [meters]')



    % Loop through all the parameters
    for ith_parameter_set = 1:N_conversions

        % Get the current parameter set
        current_parameters_type = parameters_to_convert_type_strings{ith_parameter_set};
        current_parameters = parameters_to_convert{ith_parameter_set};
       
        % Plot the inputs
        fcn_geometry_plotGeometry(current_parameters_type,current_parameters);
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

    temp_axis = axis;

    subplot(1,2,2);
    hold on;
    grid on;
    axis equal;
    title('Revised');
    xlabel('X [meters]');
    ylabel('Y [meters]')

    % Plot the outputs
    % Loop through all the parameters
    for ith_parameter_set = 1:N_conversions

        % Get the current parameter set
        current_parameters_type = parameters_to_convert_type_strings{ith_parameter_set};
        current_parameters = St_parameters{ith_parameter_set};
       
        % Plot the inputs
        fcn_geometry_plotGeometry(current_parameters_type,current_parameters);
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

