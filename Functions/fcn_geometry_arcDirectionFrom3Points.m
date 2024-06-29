function is_counterClockwise = fcn_geometry_arcDirectionFrom3Points(points1, points2, points3,varargin)
% fcn_geometry_arcDirectionFrom3Points calculates whether or not 3 points
% make an arc in the clockwise or counter-clockwise diretion. The direction
% is calculated by performing a cross product between the vectors from
% points1 to points3, versus the vectors from points1 to points2. The
% funtion returns the sign of the result.
%
% FORMAT:
%
% is_counterClockwise = fcn_geometry_arcDirectionFrom3Points(point1, point2, point3,(fig_num))
%
% INPUTS:
%
%      point1, point2, point3: a Nx2 vectors of point pairings. Note: this
%      function is vectorized so that multiple points can be entered
%      simultaneously
%
%      (OPTIONAL INPUTS)
%
%      fig_num: a figure number to plot results.
%
% OUTPUTS:
%
%      is_counterClockwise: an [Nx1] vector containing 1 if the
%      corresponding row of points creates a counter-clockwise arc
%      (positive) or -1 if it would create a clockwise arc (e.g. negative).
%      It returns 0 if the direction is undefined (if points are repeated,
%      colinear, exactly on start or end, etc.). To be defined, the
%      cross-product between the vector12 and vector 13 must have magnitude
%      greater than 1E-10.
%
% DEPENDENCIES:
%
%      fcn_DebugTools_checkInputsToFunctions
%      fcn_geometry_fillColorFromNumberOrName
%
% EXAMPLES:   
%
% See the script: script_test_fcn_geometry_arcDirectionFrom3Points
% for a full test suite.
%
% This function was written on 2023_12_19 by S. Brennan
% Questions or comments? sbrennan@psu.edu 

% Revision history:
% 2023_12_19 - sbrennan@psu.edu
% -- original write of the code
% 2024_01_03 - S. Brennan
% -- added fast mode option
% -- added environmental variable options
% 2024_01_08 - S. Brennan
% -- changed plotting to plot each case separately
% -- fixed bug with cross function call to force it to cross column-wise
% 2024_04_14 - S. Brennan
% -- added fcn_geometry_fillColorFromNumberOrName to plotting
% 2024_05_15 - S. Brennan
% -- fixed bug where gives wrong answer if test points are on start or end
% points, or VERY close to zero. Added check on numerical tolerance to fix
% this, returning 0 as result if points are at end.


%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==4 && isequal(varargin{end},-1))
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

if (0==flag_max_speed)
    if flag_check_inputs
        % Are there the right number of inputs?
        narginchk(3,4);

        % Check the points1 input
        fcn_DebugTools_checkInputsToFunctions(...
            points1, '2column_of_numbers');

        N_points = length(points1(:,1));

        % Check the points2 input
        fcn_DebugTools_checkInputsToFunctions(...
            points2, '2column_of_numbers',[N_points N_points]);

        % Check the points2 input
        fcn_DebugTools_checkInputsToFunctions(...
            points2, '2column_of_numbers',[N_points N_points]);
    end
end

N_points = length(points1(:,1));

% Does user want to show the plots?
flag_do_plots = 0;
if (4 == nargin) && (0==flag_max_speed)
    temp = varargin{end};
    if ~isempty(temp)
        fig_num = temp;
        flag_do_plots = 1;
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
 
% The code below is set up to be vectorized. The method is to create a
% vector from points1 to points2 and to points3. If the cross product is
% positive, then the arc is CCW.

projection_points1_to_points2 = points2 - points1;
projection_points1_to_points3 = points3 - points1;

% Do the cross product from segment 1-2 to segment 1-3
cross_product = cross([projection_points1_to_points2 zeros(N_points,1)],[projection_points1_to_points3 zeros(N_points,1)],2);
cross_product_direction = cross_product(:,3);

% If barely positive or negative, within numerical tolerance of 1E-10, then
% it is indeterminate and force it to be zero
cross_product_direction(abs(cross_product_direction)<1E-10) = 0;

is_counterClockwise = sign(cross_product_direction);

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
    
    %Prep the figure
    temp_h = figure(fig_num);
    flag_rescale_axis = 0;
    if isempty(get(temp_h,'Children'))
        flag_rescale_axis = 1;
    end        

    tiledlayout('flow')

    % Calculate the location for text
    half_point = (points1+points3)/2;
    half_toward_point2 = (half_point+points2)/2;
    text_locations = half_toward_point2; % +half_ortho_distance.*unit_orthogonal_points1_to_points3;
    difference = sum((half_point-points2).^2,2).^0.5;

    % plot all the results    
    for i_fit = 1:N_points    

        nexttile;

        hold on;
        grid on;
        grid minor;
        axis equal;

        if is_counterClockwise(i_fit)==1
            label_clockwise_or_counterclockwise = 'Counter-clockwise';
        else
            label_clockwise_or_counterclockwise = 'Clockwise';
        end

        title(sprintf('Point set %.0d: %s', i_fit, label_clockwise_or_counterclockwise),'Interpreter','none');


        % Plot the inputs
        plot(points1(i_fit,1),points1(i_fit,2),'g.','MarkerSize',30);  % Plot points1
        plot(points2(i_fit,1),points2(i_fit,2),'b.','MarkerSize',20);  % Plot points2
        plot(points3(i_fit,1),points3(i_fit,2),'r.','MarkerSize',10);  % Plot points3


        current_color = fcn_geometry_fillColorFromNumberOrName(i_fit);


        % Plot the connecting lines
        plot([points1(i_fit,1) points3(i_fit,1)],[points1(i_fit,2) points3(i_fit,2)],'-','LineWidth',3,'Color',current_color);

        % Plot the midpoint lines
        plot([half_point(i_fit,1) points2(i_fit,1)],[half_point(i_fit,2) points2(i_fit,2)],'-','LineWidth',3,'Color',current_color);

        % Label the result
        nudge = 0.1*difference(i_fit);
        text(text_locations(i_fit,1)+nudge,text_locations(i_fit,2)+nudge,sprintf('%.0d',is_counterClockwise(i_fit)));

        % Make axis slightly larger?
        if flag_rescale_axis
            temp = axis;
            %     temp = [min(points(:,1)) max(points(:,1)) min(points(:,2)) max(points(:,2))];
            axis_range_x = temp(2)-temp(1);
            axis_range_y = temp(4)-temp(3);
            percent_larger = 0.3;
            axis([temp(1)-percent_larger*axis_range_x, temp(2)+percent_larger*axis_range_x,  temp(3)-percent_larger*axis_range_y, temp(4)+percent_larger*axis_range_y]);
        end
    end


end

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

