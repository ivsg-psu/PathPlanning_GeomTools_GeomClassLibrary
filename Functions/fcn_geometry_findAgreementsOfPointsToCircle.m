function agreement_indicies = fcn_geometry_findAgreementsOfPointsToCircle(points, circleCenter, circleRadius, transverse_tolerance, varargin) 
% fcn_geometry_findAgreementsOfPointsToCircle
% Given a set of XY points, a circleCenter, and a circleRadius, finds the indicies of the points that are within a
% transverse_tolerance distance away from the circle. 
%
% Format:
% agreement_indicies = fcn_geometry_findAgreementsOfPointsToCircle(points, circleCenter, circleRadius, transverse_tolerance, (fig_num))
%
% INPUTS:
%
%      points: a Nx2 vector where N is the number of points, but at least 2
%      rows.
%
%      circleCenter: a 1x2 vector that is a the [x y] coordinates of the
%      circle center to test.
%
%      circleRadius: the radius of the circle to test
%
%      transverse_tolerance: the orthogonal distance between the points and
%      the linear vector fit that indicate whether a point "belongs" to the
%      fit. A point belongs to the fit if the transverse distance is less
%      than or equal to the transverse_tolerance
%
%      (OPTIONAL INPUTS)
%
%      fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed.
%
% OUTPUTS:
%
%      agreement_indicies: the indicies of the points that are within
%      agreement of the best-fit parameters, given the transverse and
%      station tolerance settings. If a station_tolerance is specified, the
%      indicies are sorted in station-increasing order.
% 
% DEPENDENCIES:
%      fcn_DebugTools_checkInputsToFunctions
%      fcn_geometry_findPointsInSequence
%
% EXAMPLES:
%
% See the script: script_test_fcn_geometry_findAgreementsOfPointsToCircle
% for a full test suite.

% This function was written on 2024_01_15 by S. Brennan
% Questions or comments? sbrennan@psu.edu

% Revision history:
% 2023_01_15/16 - S. Brennan
% -- wrote the code

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

        % Check the points input to be length greater than or equal to 2
        fcn_DebugTools_checkInputsToFunctions(...
            points, '2column_of_numbers',[2 3]);

        % Check the circleCenter input to be [1x2]
        fcn_DebugTools_checkInputsToFunctions(...
            circleCenter, '2column_of_numbers',[1 1]);

        % Check the circleRadius input to be [1x1] positive number
        fcn_DebugTools_checkInputsToFunctions(...
            circleRadius, 'positive_1column_of_numbers',1);

        % Check the transverse_tolerance input is a positive single number
        fcn_DebugTools_checkInputsToFunctions(transverse_tolerance, 'positive_1column_of_numbers',1);

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

flag_use_domain_method = 0;
% Find agreement domain for a circle fit independent of station_tolerance
if flag_use_domain_method

    % Create a domain by doing a large range of angles across an inner and
    % outer arc that spans the test area
    angles = (0:1:359.9999)'*pi/180;
    inner_radius = max(0,(circleRadius - transverse_tolerance));
    outer_radius = circleRadius + transverse_tolerance;
    inner_arc = inner_radius*[cos(angles) sin(angles)] + ones(length(angles(:,1)),1)*circleCenter;
    outer_arc = outer_radius*[cos(angles) sin(angles)] + ones(length(angles(:,1)),1)*circleCenter;
    best_fit_domain_box = [inner_arc; flipud(outer_arc)];

    % Find the points within the domain using isinterior
    domainPolyShape = polyshape(best_fit_domain_box(:,1),best_fit_domain_box(:,2),'Simplify',false,'KeepCollinearPoints',true);
    indicies_in_transverse_agreement = find(isinterior(domainPolyShape,points) == 1);

else
    % Distance of all the input points from the center
    point_radii = sum((points - circleCenter).^2,2).^0.5;

    % Absolute Error to find the indices in agreement
    absolute_radial_error = abs(point_radii - circleRadius);

    % Indices in transverse agreement
    indicies_in_transverse_agreement_binary_form = (absolute_radial_error <= transverse_tolerance)' ;
    indicies_in_transverse_agreement = find(indicies_in_transverse_agreement_binary_form==1);
end

agreement_indicies = indicies_in_transverse_agreement;



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
    axis equal;

    title('Points and agreements, plotted in point-space');
    xlabel('X [meters]');
    ylabel('Y [meters]')


    % Plot the input points
    plot(points(:,1),points(:,2),'k.','MarkerSize',30)

    % Plot the boundaries of the circle fit
    fcn_geometry_plotCircle(circleCenter, circleRadius, 'b-',fig_num); 
    fcn_geometry_plotCircle(circleCenter, circleRadius-transverse_tolerance, 'r-',fig_num); 
    fcn_geometry_plotCircle(circleCenter, circleRadius+transverse_tolerance, 'r-',fig_num); 
    plot(circleCenter(:,1),circleCenter(:,2),'b+','MarkerSize',15);


    % Plot the radial agreement points
    plot(points(agreement_indicies,1),points(agreement_indicies,2),'r.','MarkerSize',15);

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



