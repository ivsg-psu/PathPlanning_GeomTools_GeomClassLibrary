function [revised_fitSequence_types, revised_fitSequence_parameters, max_feasibility_distance] = ...
    fcn_geometry_alignGeometriesInSequence(input_types, input_parameters, threshold, varargin)
%% fcn_geometry_alignGeometriesInSequence
% Given the results of regression fits that are used to fit a set of data
% in sequence, this function proceeds from the first fit to the second,
% then second to third, etc. and matches the fits together to achieve
% levels of continuity between fits.
%
% The following matches are supported with either C0, C1, and C2
% continuity, and user-defined tolerances:
%
% * connections from line segments to arcs using
%   fcn_geometry_alignSegmentArc 
%
% * connections from arcs to line segments using
%   fcn_geometry_alignArcSegment
%
% * connections from line arcs to arcs, using
%   fcn_geometry_alignArcArc 
%
% Format:
% revised_fitSequence_parameters = ...
% fcn_geometry_alignGeometriesInSequence(fitSequence_bestFitType, fitSequence_parameters, threshold, (continuity_level), (fig_num))
%
% INPUTS:
%      fitSequence_bestFitType: a cell array of length N that contains
%      identifier strings labeling the fit type of each fit in the
%      sequence.
%
%      fitSequence_parameters: a cell array of length N that contains a
%      vector of parameters used for each fit representation
%
%      threshold: the user-defined threshold representing the allowable
%      error in sequence match. If the change in the fitted arc's position
%      is less than the allowable error, then the fit is permitted and
%      proceeds. If the change is larger than the allowable error, the fit
%      process is stopped.
%
%      (OPTIONAL INPUTS)
%
%      continuity_level: the level of continuity requested between the
%      geometries, either 0, 1, or 2. The default is 2. 
% 
%           The levels of alignment continuity, from 0 to 2, mean the
%           following:
% 
%           C0 continuous: the ith+1 segment is forced to intersect the ith
%           segment
% 
%           C1 continous: the ith+1 segment is C0 continous with the ith
%           segment, and the tangent vector of the ith+1 segment is forced
%           to match the tangent vector of the ith segment at the
%           intersection point. In other words, the derivative dt/ds in
%           St-coordinates are matched. Or, in XY coordinates, the radius
%           of ith fit at point of contact is matched to the radius of the
%           ith+1 fit at the same contact point. NOTE: a function can be C1
%           continuous in St coordinates but NOT C1 continuous in XY
%           coordinates. For example, an arc joining a vertical line has an
%           undefined slope at the vertical line, but still be tangent to
%           the circle.
% 
%           C2 continous: the ith+1 segment is C0 and C1 continous with the
%           ith segment, and the curvature (1/R) of the ith+1 segment is
%           forced to match the curvature of the ith segment at the
%           intersection point. In other words, the 2nd derivative dt/ds in
%           St-coordinates are matched such that there are no
%           discontinuities in the curvature. For vehicles, C2 continuity
%           is usually required for high-speed roads because the steering
%           angle is kinematically (but not dynamically) associated with
%           the radius of the local path. Enforcing C2 curvature thus
%           forces the curve fit to produce a path such that the steering
%           angle of a vehicle does not have to instantaneously change
%           from one position to another in order for a vehicle to remain
%           centered on a path.
%
%      fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed.
%
% OUTPUTS:
%
%      revised_fitSequence_types: a cell array of length N that contains
%      identifier strings labeling the fit type of each fit in the
%      sequence.
%
%      revised_fitSequence_parameters: the parameters for each of the N
%      fits such that they are aligned.
%
%      max_feasibility_distance: the minimum tolerance that would cause
%      feasibilty in fitting.
%
% DEPENDENCIES:
%
%      fcn_geometry_fitArcRegressionFromHoughFit
%
% EXAMPLES:
%
% See the script: script_test_fcn_geometry_alignGeometriesInSequence
% for a full test suite.
%
% This function was written on 2024_04_19 by S. Brennan
% Questions or comments? sbrennan@psu.edu

% Revision history:
% 2024_04_19 - S Brennan
% -- wrote the code
% 2024_06_21 - Sean Brennan
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
% -- changed line parameter format to new standard:
%             [
%              base_point_x, 
%              base_point_y, 
%              heading,
%             ]
% -- changed segment parameter format to new standard:
%             [
%              base_point_x, 
%              base_point_y, 
%              heading,
%              s_Length,
%             ]
% 2024_06_21 - Sean Brennan
% -- added continuity_level as an input option


%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==5 && isequal(varargin{end},-1))
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
        narginchk(3,5);

    end
end

% Does user want to specify continuity_level?
continuity_level = 2;
if (4<=nargin)
    temp = varargin{1};
    if ~isempty(temp)
        continuity_level = temp;
    end
end


% Does user want to specify fig_num?
fig_num = []; % Default is to have no figure
flag_do_plots = 0;
if  (0==flag_max_speed) && (5<= nargin)
    temp = varargin{end};
    if ~isempty(temp)
        fig_num = temp;
        flag_do_plots = 1;

        % Does user want to specify animation_figure_handles?
        % flag_plot_subfigs = 0;

        if length(temp)>1
            fig_num           = temp(1);
            % h_plotPoints      = temp(2);
            % h_plotPercentage  = temp(3);
            % h_plotFitShape    = temp(4);
            % flag_plot_subfigs = 1;
        end
    end
end


%% Solve for the circle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _
%  |  \/  |     (_)
%  | \  / | __ _ _ _ __
%  | |\/| |/ _` | | '_ \
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if flag_do_plots
        % Plot the results in the given figure number
        temp_h = figure(fig_num);
        flag_rescale_axis = 0;
        if isempty(get(temp_h,'Children'))
            flag_rescale_axis = 1;
        end

        % Plot the input geometries 
        hold on;
        grid on;
        axis equal;
        xlabel('X [meters]');
        ylabel('Y [meters]');

        segment_length = [];
        format_string = sprintf(' ''-'',''Color'',[0.6 0.6 0.6],''LineWidth'',7 ');
        fcn_geometry_plotFitSequences(input_types, input_parameters, segment_length, format_string, (fig_num));


        % Make axis slightly larger?
        if flag_rescale_axis
            temp = axis;
            %     temp = [min(points(:,1)) max(points(:,1)) min(points(:,2)) max(points(:,2))];
            axis_range_x = temp(2)-temp(1);
            axis_range_y = temp(4)-temp(3);
            percent_larger = 0.3;
            axis([temp(1)-percent_larger*axis_range_x, temp(2)+percent_larger*axis_range_x,  temp(3)-percent_larger*axis_range_y, temp(4)+percent_larger*axis_range_y]);
        end

        good_axis = axis;

end

fitSequence_bestFitType = input_types;
fitSequence_parameters = input_parameters;

% How many fits are we aligning?
NfitsInSequence = length(fitSequence_bestFitType);

% Fill in starter values for the output parameters by using an empty cell
% array
N_revisedFits                     = 0;
revised_fitSequence_types{1}      = '';
revised_fitSequence_parameters{1} = [];
max_feasibility_distance          = -inf;
lastFit_type                      = fitSequence_bestFitType{1};
lastFit_parameters                = fitSequence_parameters{1}; 

in_boundary_margin = 0.01;

% Loop through fits, connecting them together
for ith_fit = 1:NfitsInSequence-1

    if ith_fit==5
        disp('Stop here');
    end

    fprintf(1,'Performing fit %.0d of %.0d\n',ith_fit,NfitsInSequence-1);

    current_fit_type = fcn_INTERNAL_covertComplexShapeNamesToSimpleNames(lastFit_type);
    next_fit_type    = fcn_INTERNAL_covertComplexShapeNamesToSimpleNames(fitSequence_bestFitType{ith_fit+1});

    current_fit_parameters = lastFit_parameters;
    next_fit_parameters    = fitSequence_parameters{ith_fit+1};

    switch current_fit_type
        case 'segment'
            [revised_subSequence_types, revised_sequence_parameters, flag_fit_failed, feasibility_distance] = fcn_INTERNAL_alignSegmentToX(current_fit_parameters, next_fit_parameters, next_fit_type, continuity_level, threshold, in_boundary_margin);            
        case 'arc'
            [revised_subSequence_types, revised_sequence_parameters, flag_fit_failed, feasibility_distance] = fcn_INTERNAL_alignArcToX(current_fit_parameters, next_fit_parameters, next_fit_type, continuity_level, threshold, in_boundary_margin);
        otherwise
            warning('on','backtrace');
            warning('An error will be thrown at this point due to missing code.');
            error('Alignments are not yet supported for curves from fit type: %s',current_fit_type);
    end

    max_feasibility_distance = max(max_feasibility_distance,feasibility_distance);
    
    if flag_fit_failed
        break;
    end

    % Save the results into the cumulative cell array, revised_fitSequence,
    % that is storing all the fits. 
    N_subSequences = length(revised_subSequence_types);
    for ith_result = 1:N_subSequences-1
        N_revisedFits = N_revisedFits+1;
        revised_fitSequence_types{N_revisedFits}      = revised_subSequence_types{ith_result}; %#ok<AGROW>
        revised_fitSequence_parameters{N_revisedFits} = revised_sequence_parameters{ith_result}; %#ok<AGROW>
    end

    % Save the last fit - this is used at the start of the next loop
    % through
    lastFit_type       = revised_subSequence_types{N_subSequences};
    lastFit_parameters = revised_sequence_parameters{N_subSequences};


    if flag_do_debug
        
        debug_fig_num_iterated = ith_fit+56575;
        figure(debug_fig_num_iterated);

        sgtitle('Debugging')
        subplot(1,2,2);
        cla;
        axis equal;
        fcn_geometry_plotFitSequences(revised_subSequence_types, revised_sequence_parameters,(debug_fig_num_iterated));
        temp_axis = axis;
        title('Results from last joint alignment');

        subplot(1,2,1);
        cla;
        axis equal
        fcn_geometry_plotFitSequences(input_types, input_parameters,(debug_fig_num_iterated));        
        fcn_geometry_plotFitSequences(revised_fitSequence_types, revised_fitSequence_parameters,(debug_fig_num_iterated));
        axis(temp_axis);
        title('Cumulative alignments')

    elseif flag_do_plots
        figure(fig_num);

        fcn_geometry_plotFitSequences(revised_fitSequence_types, revised_fitSequence_parameters,[],[],(fig_num));

        % Match axis
        axis(good_axis);
    end

end % Ends looping through fits

% Save the last result
N_revisedFits = N_revisedFits+1;
revised_fitSequence_types{N_revisedFits}      = lastFit_type; 
revised_fitSequence_parameters{N_revisedFits} = lastFit_parameters; 

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

    % if flag_plot_subfigs
    %
    %     figure(get(h_plotFitShape.Parent.Parent, 'Number'));
    %
    %     % Plot the results in the subplot
    %     flag_rescale_axis = 0;
    %
    %     % Match subplot 4 axis with that from subplot 1
    %     subplot(2,2,1);
    %     original_axis = axis;
    %
    % else
    % Plot the results in the given figure number
    temp_h = figure(fig_num);
    flag_rescale_axis = 0;
    if isempty(get(temp_h,'Children'))
        flag_rescale_axis = 1;
    end
    % end


    % Plot the input geometries in grey
    hold on;
    grid on;
    axis equal;
    xlabel('X [meters]');
    ylabel('Y [meters]');

    segment_length = [];
    format_string = sprintf(' ''-'',''Color'',[0.6 0.6 0.6],''LineWidth'',7 ');
    fcn_geometry_plotFitSequences(input_types, input_parameters, segment_length, format_string, (fig_num));


    % Make axis slightly larger?
    if flag_rescale_axis
        temp = axis;
        %     temp = [min(points(:,1)) max(points(:,1)) min(points(:,2)) max(points(:,2))];
        axis_range_x = temp(2)-temp(1);
        axis_range_y = temp(4)-temp(3);
        percent_larger = 0.3;
        axis([temp(1)-percent_larger*axis_range_x, temp(2)+percent_larger*axis_range_x,  temp(3)-percent_larger*axis_range_y, temp(4)+percent_larger*axis_range_y]);
    end

    good_axis = axis;

    % Plot the aligned geometries on the left
    fcn_geometry_plotFitSequences(revised_fitSequence_types, revised_fitSequence_parameters,[],[],(fig_num));

    % Match axis
    axis(good_axis);

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

%% fcn_INTERNAL_covertComplexShapeNamesToSimpleNames
function simple_name_string = fcn_INTERNAL_covertComplexShapeNamesToSimpleNames(complex_name_string)

switch lower(complex_name_string)
    case {'arc','line','segment','spiral','','none','circle'}
        simple_name_string = complex_name_string;
    case {'regression arc'}
        simple_name_string = 'arc';
    case {'vector regression segment fit'}
        simple_name_string = 'segment';
    otherwise
        warning('on','backtrace');
        warning('An error will be thrown due to unrecognized fitting name type, inside fcn_INTERNAL_covertComplexShapeNamesToSimpleNames.');
        error('Unrecognized fit string: %s', complex_name_string);
end

end % Ends fcn_INTERNAL_covertComplexShapeNamesToSimpleNames

%% fcn_INTERNAL_alignArcToX
function [revised_sequence_types, revised_sequence_parameters, flag_fit_failed, feasibility_distance] = fcn_INTERNAL_alignArcToX(arc_parameters, X_parameters, X_fitType, continuity_level, threshold, in_boundary_margin)

if length(threshold(1,:))>1
    transverse_threshold = threshold(1,2);
else
    transverse_threshold = threshold(1,1);
end

% Initialize values
feasibility_distance = 0;

switch X_fitType
    case 'segment'
        % Arc to Segment
        
        segment_parameters = X_parameters;

        if 2==continuity_level
            %%%%
            % Check feasibility
            % Format:
            % [flag_is_feasible, feasibility_distance, closest_feasible_line_parameters] = ...
            % fcn_geometry_isC2FeasibleArcToLine( circle_parameters, line_parameters, (threshold), (in_boundary_margin), (fig_num));

            % Is a C2 solution feasible?
            [flag_is_feasible, feasibility_distance, closest_feasible_segment_parameters] = ...
                fcn_geometry_isC2FeasibleArcToLine( arc_parameters, segment_parameters, (transverse_threshold), (in_boundary_margin), (-1));

            % If not C2, then is a C1 solution feasible?
            if 0==flag_is_feasible
                [flag_is_feasible, feasibility_distance, closest_feasible_segment_parameters] = ...
                    fcn_geometry_isC2FeasibleArcToLine( arc_parameters, segment_parameters, (transverse_threshold), (in_boundary_margin), (-1));
                if 0==flag_is_feasible
                    error('Curves encountered that are neither C2 or C1 feasible - this should not be possible.');
                    revised_arc_parameters = nan*arc_parameters;
                end
                continuity_level = 1;

            end
        else
            closest_feasible_segment_parameters = segment_parameters;
            flag_is_feasible = 1;
        end

        % If possible, try to calculate the revised parameters
        if 1==flag_is_feasible
            %%%%
            % Calculate parameters
            [revised_arc_parameters, revised_X_parameters, revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = fcn_geometry_alignArcSegment(...
                arc_parameters, closest_feasible_segment_parameters, (threshold), (continuity_level), (-1));
        end

    case 'arc'
        % Arc to Arc2

        arc2_parameters = X_parameters;
        
        if 2==continuity_level
            %%%%
            % Check feasibility
            % Format:
            % [flag_is_feasible, feasibility_distance, closest_feasible_arc2_parameters] = ...
            % fcn_geometry_isC2FeasibleArcToArc(arc_parameters, arc2_parameters, (threshold), (in_boundary_margin), (fig_num));

            % Is a C2 solution feasible?
            [flag_is_feasible, feasibility_distance, closest_feasible_arc2_parameters] = ...
            fcn_geometry_isC2FeasibleArcToArc(arc_parameters, arc2_parameters, (transverse_threshold), (in_boundary_margin), (-1));

            % If not C2, then is a C1 solution feasible?
            if 0==flag_is_feasible
                [flag_is_feasible, feasibility_distance, closest_feasible_arc2_parameters] = ...
                    fcn_geometry_isC1FeasibleArcToArc(arc_parameters, arc2_parameters, (transverse_threshold), (in_boundary_margin), (-1));
                if 0==flag_is_feasible
                    error('Curves encountered that are neither C2 or C1 feasible - this should not be possible.');
                    revised_arc_parameters = nan*arc_parameters;
                end
                continuity_level = 1;
            end
            
        else
            closest_feasible_arc2_parameters = arc2_parameters;
            flag_is_feasible = 1;
        end

        % If possible, try to calculate the revised parameters
        if 1==flag_is_feasible
            %%%%
            % Calculate parameters
            [revised_arc_parameters, revised_X_parameters, revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = fcn_geometry_alignArcArc(...
                arc_parameters, closest_feasible_arc2_parameters, (threshold), (continuity_level), (2689));
        end

    otherwise
        warning('on','backtrace');
        warning('An error will be thrown at this point due to missing code.');
        error('Alignments are not yet supported for curves from fit type: %s',current_fit_type);
end

% Check to see if any of the functions return NaN values. If so, this is a
% failed fit
flag_fit_failed = 0;
if any(isnan(revised_arc_parameters),'all')
    flag_fit_failed = 1;
elseif any(isnan(revised_X_parameters),'all')
    flag_fit_failed = 1;
else
    if ~any(isnan(revised_intermediate_geometry_join_parameters),'all')
        revised_sequence_types{1}      = 'arc';
        revised_sequence_types{2}      = revised_intermediate_geometry_join_type;
        revised_sequence_types{3}      = X_fitType;

        revised_sequence_parameters{1} = revised_arc_parameters;
        revised_sequence_parameters{2} = revised_intermediate_geometry_join_parameters;
        revised_sequence_parameters{3} = revised_X_parameters;
    else
        revised_sequence_types{1}      = 'arc';
        revised_sequence_types{2}      = X_fitType;

        revised_sequence_parameters{1} = revised_arc_parameters;
        revised_sequence_parameters{2} = revised_X_parameters;

    end
end

if 1==flag_fit_failed
    revised_sequence_types{1} = '';
    revised_sequence_parameters{1} = nan;
end
end % Ends fcn_INTERNAL_alignArcToX


%% fcn_INTERNAL_alignSegmentToX
function [revised_sequence_types, revised_sequence_parameters, flag_fit_failed, feasibility_distance] = fcn_INTERNAL_alignSegmentToX(segment_parameters, X_parameters, X_fitType, continuity_level, threshold, in_boundary_margin)

if length(threshold(1,:))>1
    transverse_threshold = threshold(1,2);
else
    transverse_threshold = threshold(1,1);
end

% Initialize values
feasibility_distance = 0;

switch X_fitType
    case 'segment'
        % Segment to Segment
        warning('on','backtrace');
        warning('An error will be thrown at this point due to missing code.');
        error('Alignments from segment to segment are not yet supported.');

    case 'arc'
        % Segment to Arc
        
        arc_parameters = X_parameters;

        if 2==continuity_level
            %%%%
            % Check feasibility
            % Format:
            % [flag_is_feasible, feasibility_distance, closest_feasible_circle_parameters] = ...
            % fcn_geometry_isC2FeasibleLineToArc(line_parameters, circle_parameters, (threshold), (in_boundary_margin), (fig_num));
            [flag_is_feasible, feasibility_distance, closest_feasible_arc_parameters] = ...
                fcn_geometry_isC2FeasibleLineToArc(segment_parameters, arc_parameters, (transverse_threshold), (in_boundary_margin), (-1));

            % Is a solution feasible?
            if 0==flag_is_feasible
                revised_segment_parameters = nan*segment_parameters;
            end
        else
            closest_feasible_arc_parameters = arc_parameters;
            flag_is_feasible = 1;
        end

        % If possible, try to calculate the revised parameters
        if 1==flag_is_feasible
            %%%%
            % Calculate parameters
            % Format:
            % [revised_segment_parameters, revised_arc_parameters, revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters]  = ...
            % fcn_geometry_alignSegmentArc(segment_parameters, arc_parameters, (threshold), (continuity_level),  (fig_num))
            %
            % Call function
            [revised_segment_parameters, revised_X_parameters, revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = fcn_geometry_alignSegmentArc(...
                segment_parameters, closest_feasible_arc_parameters, (threshold), (continuity_level), (-1));
        end

    otherwise
        warning('on','backtrace');
        warning('An error will be thrown at this point due to missing code.');
        error('Alignments are not yet supported for curves from fit type: %s',current_fit_type);
end

% Check to see if any of the functions return NaN values. If so, this is a
% failed fit
flag_fit_failed = 0;
if any(isnan(revised_segment_parameters),'all')
    flag_fit_failed = 1;
elseif any(isnan(revised_X_parameters),'all')
    flag_fit_failed = 1;
else
    if ~any(isnan(revised_intermediate_geometry_join_parameters),'all')
        revised_sequence_types{1}      = 'segment';
        revised_sequence_types{2}      = revised_intermediate_geometry_join_type;
        revised_sequence_types{3}      = X_fitType;

        revised_sequence_parameters{1} = revised_segment_parameters;
        revised_sequence_parameters{2} = revised_intermediate_geometry_join_parameters;
        revised_sequence_parameters{3} = revised_X_parameters;
    else
        revised_sequence_types{1}      = 'segment';
        revised_sequence_types{2}      = X_fitType;

        revised_sequence_parameters{1} = revised_segment_parameters;
        revised_sequence_parameters{2} = revised_X_parameters;

    end
end
end % Ends fcn_INTERNAL_alignSegmentToX



