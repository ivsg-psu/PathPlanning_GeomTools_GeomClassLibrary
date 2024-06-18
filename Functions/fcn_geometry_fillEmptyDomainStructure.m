function emptyDomain = fcn_geometry_fillEmptyDomainStructure(varargin)
% fcn_geometry_fillEmptyDomainStructure
% Fills in an empty domain-type structure with the minimum fields. These
% fields include:
%
% Format:
% emptyDomain = fcn_geometry_fillEmptyDomainStructure((fig_num))
%
% INPUTS:
%
%      (OPTIONAL INPUTS)
%
%      fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed.
%
% OUTPUTS:
%
%      emptyDomain: a domain-type structure containing the following
%      fields:
%
%      emptyDomain.best_fit_type = 'empty';
%         
%            This is a string used to denote the fit type: 
% 
%            'empty' - no fit yet done
%
%            'unfitted' - points left over after a fit is done 
%            
%            'Hough line' - A Hough fit of line type
%            
%            'Hough segment' - A Hough fit of line segment type
%            
%            'Hough circle' - A Hough fit of a circle type
%            
%            'Hough arc' - A Hough fit of arc type
%
%            'Vector regression line fit' - a vector fit of line data,
%            using regression
%            
%            'Vector regression segment fit' - a vector fit of a line
%            segment, using regression
%
%            'spiral' - an Euler spiral curve such that the curvature is
%            linearly changing versus the station coordinate.
%    
%            'Hough cubic polynomial' - A Hough fit of cubic polynomial
%            type
%            
%      emptyDomain.points_in_domain = [nan nan];
%
%            These are the points in [N x 2] or [N x 3] format that were
%            used to create the domain, particularly the statistics of the
%            domain
%   
%      emptyDomain.best_fit_parameters = nan;
% 
%            These are the parameters of the fit. They are defined
%            in each fit type in the following manner:
%
%            'empty' - no fit yet done, so [nan]
%
%            'unfitted' - points left over after a fit is done, so [nan]
%            
%            'Hough line' - 
%
%             [unit_projection_vector_x,
%              unit_projection_vector_y,
%              base_point_x, 
%              base_point_y, 
%             ]
%            
%            'Hough segment' - 
%
%             [unit_projection_vector_x,
%              unit_projection_vector_y,
%              base_point_x, 
%              base_point_y, 
%              station_distance_min,
%              station_distance_max,
%             ]
%
%            'Hough circle' - 
%               [circleCenter_x.
%                circleCenter_y,
%                radius]
%            
%            'Hough arc' - 
%
%               [circleCenter_x.
%                circleCenter_y,
%                radius,
%                start_angle_in_radians, 
%                end_angle_in_radians,
%                flag_this_is_a_circle
%                flag_arc_is_counterclockwise
%               ] 
%
%            'Vector regression line fit' - 
%
%             [unit_projection_vector_x,
%              unit_projection_vector_y,
%              base_point_x, 
%              base_point_y, 
%             ]
%
%            'Vector regression segment fit' -

%             [unit_projection_vector_x,
%              unit_projection_vector_y,
%              base_point_x, 
%              base_point_y, 
%              station_distance_min,
%              station_distance_max,
%             ]
%
%            'Regression circle' - 
%               [circleCenter_x.
%                circleCenter_y,
%                radius]
%            
%            'Regression arc' - 
%
%               [circleCenter_x.
%                circleCenter_y,
%                radius,
%                start_angle_in_radians, 
%                end_angle_in_radians,
%                flag_this_is_a_circle
%                flag_arc_is_counterclockwise
%               ] 
%            
%            'spiral' - 
%
%               [spiralLength,  % the s-coordinate length allowed
%                h0,  % The initial heading
%                x0,  % The initial x value
%                y0,  % The initial y value
%                K0,  % The initial curvature
%                Kf   % The final curvature
%              ] 
%
%            'Hough Cubic Polynomial' - 
% 
%             [ coeff_x^3, % cubic term coeffiecient
%               coeff_x^2, % square term coeffiecient
%               coeff_x^1, % linear term coeffiecient
%               coeff_x^0, % constant term coeffiecient
%              ]
%
%           
%% FIXED on 2024_06_14
%            These are the parameters of the fit. They are defined
%            in each fit type in the following manner:
%
%            'empty' - no fit yet done, so [nan]
%
%            'unfitted' - points left over after a fit is done, so [nan]
%            
%            'point' -
%             [
%              base_point_x, 
%              base_point_y, 
%             ]
%
%            'Hough line', 'Vector regression line fit' - 
%
%             [
%              base_point_x, 
%              base_point_y, 
%              heading,
%             ]
%            
%            'Hough segment', 'Vector regression segment fit' - 
%
%             [
%              base_point_x, 
%              base_point_y, 
%              heading,
%              s_Length,
%             ]
%
% DONE
%            'Hough circle', 'Regression circle' - 
%               [circleCenter_x,
%                circleCenter_y,
%                radius]
% 
% DONE
%            'Hough arc', 'Regression arc' - 
%
%               [circleCenter_x.
%                circleCenter_y,
%                radius,
%                start_angle_in_radians, 
%                end_angle_in_radians,
%                flag_arc_is_counterclockwise
%                flag_this_is_a_circle
%               ] 
%
% 
% DONE
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
%
% DONE
%            'Hough Cubic Polynomial', 'Polyfit cubic polynomial' - 
% 
%             [ 
%              base_point_x, 
%              base_point_y, 
%              flag_is_forward_x,
%              x_Length,
%              coeff_x^0, % cubic term coeffiecient
%              coeff_x^1, % square term coeffiecient
%              coeff_x^2, % linear term coeffiecient
%              coeff_x^3, % constant term coeffiecient
%             ]
%
%            
%      emptyDomain.best_fit_domain_box = [];
%
%            These are the points in [N x 2] format that define the outer
%            boundary of the best-fit domain, generally. For regression
%            fits, it is the 2-sigma box. For Hough fits, it is the
%            bounding box used for Hough voting.
%
%      emptyDomain.best_fit_1_sigma_box = [];
%
%            These are the points in [N x 2] format that define the outer
%            boundary of the domain which is the fit, +/- 1 sigma.
%
%      emptyDomain.best_fit_2_sigma_box = [];
%
%            These are the points in [N x 2] format that define the outer
%            boundary of the domain which is the fit, +/- 2 sigma.
%
%      emptyDomain.best_fit_3_sigma_box = [];
%
%            These are the points in [N x 2] format that define the outer
%            boundary of the domain which is the fit, +/- 3 sigma.

% DEPENDENCIES:
%      
%      (none)
%
% EXAMPLES:
%
% See the script: script_test_fcn_geometry_fillEmptyDomainStructure
% for a full test suite.

% This function was written on 2024_01_15 by S. Brennan
% Questions or comments? sbrennan@psu.edu

% Revision history:
% 2024_01_15 - S. Brennan
% -- wrote the code
% 2024_05_15 - Aneesh Batchu
% -- added comments for 'Hough cubic polynomial'

%% Debugging and Input checks
 
% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==1 && isequal(varargin{end},-1))
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
    if flag_check_inputs == 1
        % Are there the right number of inputs?
        narginchk(0,1);

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
if 1<= nargin && 0==flag_max_speed
    temp = varargin{end};
    if ~isempty(temp)
        fig_num = temp; %#ok<NASGU>
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
emptyDomain.best_fit_type = 'empty';
emptyDomain.points_in_domain = [nan nan];
emptyDomain.best_fit_parameters = nan;
emptyDomain.best_fit_domain_box = [];
emptyDomain.best_fit_1_sigma_box = [];
emptyDomain.best_fit_2_sigma_box = [];
emptyDomain.best_fit_3_sigma_box = [];

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
%     temp_h = figure(fig_num);
%     flag_rescale_axis = 0;
%     if isempty(get(temp_h,'Children'))
%         flag_rescale_axis = 1;
%     end
% 
%     hold on;
%     grid on;
%     title('Points and maximum-vote fit, plotted in point-space');
%     xlabel('X [meters]');
%     ylabel('Y [meters]')
% 
%     % Produce the sorted list, to create the Hough plot
%     [~,sorted_indicies] = sort(agreements,'ascend');
% 
%     % Save the number of permutations
%     N_permutations = length(combos_paired(:,1));
% 
%     % Plot the results in point space
% 
%     % Plot the input points
%     plot(points(:,1),points(:,2),'k.','MarkerSize',20);
% 
%     % Plot the best-fit points and fits
%     if flag_find_only_best_agreement ==1
%         plot(points(best_agreement_indicies,1),points(best_agreement_indicies,2),'r.','MarkerSize',15);
% 
%         % Plot the fit
%         % best_fitted_parameters_unitvector_basepoint_distance_form{N_fits}
%     else
%         for ith_agreement = 1:length(best_agreement_indicies)
%             fitted_agreement_indicies = best_agreement_indicies{ith_agreement};
%             plot(points(fitted_agreement_indicies,1),points(fitted_agreement_indicies,2),'o','MarkerSize',(15+5*ith_agreement));
%         end
%     end
% 
%     % Plot all line fits NOTE: this can be confusing as these are the INPUT
%     % line segments, not the actual line segments at the end of the fit.
%     if 1==0
%         start_points = points(combos_paired(:,1),:);
%         start_angles = atan2(start_points(:,2),start_points(:,1));
%         start_radii  = rhos ./ cos(start_angles - phis);
%         start_points_calculated = start_radii.*[cos(start_angles) sin(start_angles)];
% 
%         end_points = points(combos_paired(:,2),:);
%         end_angles = atan2(end_points(:,2),end_points(:,1));
%         end_radii  = rhos ./ cos(end_angles - phis);
%         end_points_calculated = end_radii.*[cos(end_angles) sin(end_angles)];
% 
%         N_steps = 10;
%         for ith_best_line = min(N_steps,length(agreements)):-1:1
%             ith_line = sorted_indicies(end-ith_best_line+1);
%             plot_color = (ith_best_line - 1)/N_steps*[1 1 1];
%             line_size = ceil((N_steps - ith_best_line + 1)/N_steps * 5);
% 
%             plot(...
%                 [start_points_calculated(ith_line,1) end_points_calculated(ith_line,1)],...
%                 [start_points_calculated(ith_line,2) end_points_calculated(ith_line,2)], '-','LineWidth',line_size,'Color',plot_color);
%         end
% 
% 
%     end
% 
%     % Make axis slightly larger?
%     if flag_rescale_axis
%         temp = axis;
%         %     temp = [min(points(:,1)) max(points(:,1)) min(points(:,2)) max(points(:,2))];
%         axis_range_x = temp(2)-temp(1);
%         axis_range_y = temp(4)-temp(3);
%         percent_larger = 0.3;
%         axis([temp(1)-percent_larger*axis_range_x, temp(2)+percent_larger*axis_range_x,  temp(3)-percent_larger*axis_range_y, temp(4)+percent_larger*axis_range_y]);
%     end
% 
%     % Plot the Hough space results, from least to best
%     figure(fig_num+1);
%     clf;
%     hold on;
%     grid on;
% 
%     title('Hough space plot of all fits');
%     subtitle('Points are color-weighted wherein higher darkness indicates higher votes');
% 
%     % Plot the Hough fits in batches, making the points bigger and darker
%     % in each batch to emphasize which points were the best fits.
%     N_steps = 20;
%     indicies = linspace(1,N_permutations,N_steps);
%     for ith_plot = 1:N_steps-1
%         plot_indicies_start = ceil(indicies(ith_plot));
%         plot_indicies_end = floor(indicies(ith_plot+1));
% 
%         sorted_indicies_to_plot = sorted_indicies(plot_indicies_start:plot_indicies_end);
%         plot_color = (N_steps - ith_plot + 1)/N_steps*[1 1 1];
%         marker_size = ceil(ith_plot*20/N_steps);
%         plot(phis(sorted_indicies_to_plot),rhos(sorted_indicies_to_plot),'.','MarkerSize',marker_size,'Color',plot_color);
%     end
% 
%     % Plot the best-fit phis and rhos with a red dot
%     if flag_find_only_best_agreement ==1
%             plot(best_fitted_parameters_phi_rho_form(1,1),best_fitted_parameters_phi_rho_form(1,2),'r.','MarkerSize',30);
%     else
%         for ith_agreement = 1:length(best_agreement_indicies)
%             fitted_parameters = best_fitted_parameters_phi_rho_form{ith_agreement};
%             plot(fitted_parameters(1,1),fitted_parameters(1,2),'.','MarkerSize',30);
%         end
%     end
% 
% 
%     xlabel('Phi [radians]');
%     ylabel('Rho [meters]');
% 
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


