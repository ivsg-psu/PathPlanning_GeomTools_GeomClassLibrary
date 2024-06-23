function [flag_is_feasible, feasibility_distance, closest_feasible_arc2_parameters] = fcn_geometry_isArcToArcFeasible(arc1_parameters, arc2_parameters, continuity_level, varargin)
%% fcn_geometry_isArcToArcFeasible
% Given two arc geometries, checks if the arcs can be joined with C2
% continuity. If not, it finds the closest C2 feasible arc2 given arc1.
%
% Format:
%  [flag_is_feasible, feasibility_distance, closest_feasible_arc2_parameters] = fcn_geometry_isArcToArcFeasible(...
%    arc1_parameters, arc2_parameters, ...
%    (threshold), (in_boundary_margin), (fig_num))
%
% INPUTS:
%
%      arc1_parameters: a vector of arc1 parameters consistent with arc or
%      circle type geometries. See fcn_geometry_fillEmptyDomainStructure
%      for details
%
%      arc2_parameters: a vector of arc2 parameters consistent with arc or
%      circle type geometries. See fcn_geometry_fillEmptyDomainStructure
%      for details
%
%      continuity_level: the level of continuity desired in the alignment.
%      Input values include 0 for C0 continuity, 1 for C1 continuity
%      (default), or 2 for C2 continuity. For an explanation of continuity,
%      see fcn_geometry_alignGeometriesInSequence
%
%      (OPTIONAL INPUTS)
%
%      threshold: the offset, in meters, between the arcs such that this
%      offset can be removed by shifting. If the needed feasible offset is
%      larger than this value, then the flag_is_feasible is set to 0. If
%      threshold is entered as a 2x1 or 1x2, then this specifies the
%      threshold in St coordinates, e.g. first in the station direction,
%      and then in the transverse direction. For example, an entry of [3
%      0.2] would have 3 meters threshold in the station direction, but 0.2
%      meters threshold in the transverse direction. Note that only the
%      transverse distances are used to determine feasibility. The default
%      threshold is 0
%
%      in_boundary_margin: the distance, in meters, that the offset must be
%      moved within the boundary for feasibility. This margin is added
%      because feasible solutions that exact exactly on the boundary
%      (within the tolerance) may still yield solutions that are
%      technically not feasible spiral connections. Default is 0.
%
%      fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed.
%
% OUTPUTS:
%
%      flag_is_feasible: a value of 1 if C2 continuity between arc1 and
%      arc2 is feasible for the given tolerance, 0 otherwise. 
%
%      feasibility_distance: the distance between arc1 and arc2 that gives
%      feasibile C2 continuity.
%
%      closest_feasible_arc2_parameters: the closest parameters that give
%      feasible continuity between arc2 and arc1.
%
% DEPENDENCIES:
%
%      (none)
%
% EXAMPLES:
%
% See the script: script_test_fcn_geometry_isArcToArcFeasible
% for a full test suite.
%
% This function was written on 2024_06_22 by S. Brennan
% Questions or comments? sbrennan@psu.edu

% Revision history:
% 2024_06_22 - S Brennan
% -- wrote the code

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==6 && isequal(varargin{end},-1))
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
        narginchk(3,6);

    end
end


% Does user want to specify threshold?
threshold = 0;
if (3<=nargin)
    temp = varargin{1};
    if ~isempty(temp)
        threshold = temp;
    end
end

% Does user want to specify in_boundary_margin?
in_boundary_margin = 0;
if (4<=nargin)
    temp = varargin{2};
    if ~isempty(temp)
        in_boundary_margin = temp;
    end
end


% Does user want to specify fig_num?
fig_num = []; % Default is to have no figure
flag_do_plots = 0;
if  (0==flag_max_speed) && (6<= nargin)
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


%% Solve for the closest boundaries
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _
%  |  \/  |     (_)
%  | \  / | __ _ _ _ __
%  | |\/| |/ _` | | '_ \
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch continuity_level
    case 0
        error('This is not coded yet!');        
    case 1
        error('This is not coded yet!');        
    case 2
        [flag_is_feasible, feasibility_distance, closest_feasible_arc2_parameters] = fcn_geometry_isC2FeasibleArcToArc(arc1_parameters, arc2_parameters, (threshold), (in_boundary_margin), (fig_num));
    otherwise
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
    % 
    % % Plot the results in the given figure number
    % temp_h = figure(fig_num); 
    % flag_rescale_axis = 0;
    % if isempty(get(temp_h,'Children'))
    %     flag_rescale_axis = 1;
    % end
    % % end
    % 
    % % Plot the input geometries 
    % subplot(1,2,1);
    % hold on;
    % grid on;
    % axis equal;
    % xlabel('X [meters]');
    % ylabel('Y [meters]');
    % 
    % % Plot the inputs
    % segment_length = [];
    % format1_string = sprintf(' ''-'',''Color'',[0.6 0 0],''LineWidth'',7 ');
    % 
    % if 1==flag_arc1_is_arc
    %     fcn_geometry_plotGeometry('arc', arc1_parameters,segment_length, format1_string, (fig_num));
    % else
    %     fcn_geometry_plotGeometry('circle', arc1_parameters,segment_length, format1_string, (fig_num));
    % end
    % 
    % format2_string = sprintf(' ''-'',''Color'',[1 0 0],''LineWidth'',4 ');
    % if 1==flag_arc2_is_arc
    %     fcn_geometry_plotGeometry('arc', arc2_parameters,segment_length, format2_string, (fig_num));
    % else
    %     fcn_geometry_plotGeometry('circle', arc2_parameters,segment_length, format2_string, (fig_num));
    % end
    % 
    % % Plot the corrected solution?
    % if 0==flag_is_feasible
    %     format_string = sprintf(' ''-'',''Color'',[0 1 0],''LineWidth'',4 ');
    %     if flag_check_arcs
    %         fcn_geometry_plotGeometry('arc', closest_feasible_arc2_parameters,segment_length, format_string, (fig_num));
    %     else
    %         fcn_geometry_plotGeometry('circle', closest_feasible_arc2_parameters,segment_length, format_string, (fig_num));
    %     end
    % end
    % 
    % % Make axis slightly larger?
    % if flag_rescale_axis
    %     temp = axis;
    %     %     temp = [min(points(:,1)) max(points(:,1)) min(points(:,2)) max(points(:,2))];
    %     axis_range_x = temp(2)-temp(1);
    %     axis_range_y = temp(4)-temp(3);
    %     percent_larger = 0.3;
    %     axis([temp(1)-percent_larger*axis_range_x, temp(2)+percent_larger*axis_range_x,  temp(3)-percent_larger*axis_range_y, temp(4)+percent_larger*axis_range_y]);
    % end
    % 
    % % Plot the aligned geometries on the left
    % subplot(1,2,2);
    % hold on;
    % grid on;
    % axis equal;
    % xlabel('r2/r1 [meters]');
    % ylabel('d12/r1 [meters]');
    % 
    % % Plot the cases
    % plot([0 1],[1 0],'b-','LineWidth',3);
    % plot([1 2],[0 1],'m-','LineWidth',3);
    % plot([0 1],[1 2],'c-','LineWidth',3);
    % plot(query_point(1,1)/r1,query_point(1,2)/r1,'r.','MarkerSize',30);
    % plot(corrected_parameters_r2_d12(1,1)/r1,corrected_parameters_r2_d12(1,2)/r1,'g.','MarkerSize',30);
    % 
    % legend('Case1','Case2','Case3','Query','Fixed')
    % 
    % % Make axis slightly larger?
    % if flag_rescale_axis
    %     temp = axis;
    %     %     temp = [min(points(:,1)) max(points(:,1)) min(points(:,2)) max(points(:,2))];
    %     axis_range_x = temp(2)-temp(1);
    %     axis_range_y = temp(4)-temp(3);
    %     percent_larger = 0.3;
    %     axis([temp(1)-percent_larger*axis_range_x, temp(2)+percent_larger*axis_range_x,  temp(3)-percent_larger*axis_range_y, temp(4)+percent_larger*axis_range_y]);
    % end
    % 
    % % Add a title
    % if 0==flag_is_feasible
    %     sgtitle(sprintf('NOT C2 feasible, dist: %.2f',feasibility_distance));
    % else
    %     sgtitle(sprintf('C2 feasible, dist: %.2f',feasibility_distance));
    % end

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
