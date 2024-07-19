function max_radius = fcn_geometry_maxRadiusInsideBox(box_width, box_height, varargin)
%% fcn_geometry_maxRadiusInsideBox
% Given a box defined by width and height, finds the maximum radius of an
% arc that is exactly bounded by the box, e.g an arc touching the top of
% the box while entering and exiting from the bottom corners of the box.
%
% This is useful, for example, when calculating whether an arc fit of noisy
% data is actually a line or an arc. For example: assume points are fit
% with uncertainty of +/- 2 with an arc of base width 5, and assume the
% radius of the fit is 7. This is similar to assuming an arc exists within
% a box 4 high and 5 wide. If this 4x5 box allows radii to exist that are 6
% or larger, then the fit may actually be a line, not an arc. In other
% words, if a regression fit of an arc produces a radius LARGER than that
% predicted by this function, then that radius should not be trusted - the
% data is possibly a line. A measure of the quality of the fit is the ratio
% of the fitted curvature to the curvature of this function, e.g. the
% signal to noise ratio. SNR = C_fitted/C_minimum = (1/Rfitted)/(1/R_max).
% This simplifies to SNR = max_radius/fitted_radius, wherein this function
% calculates the max_radius.
% 
% Format: 
% max_radius = fcn_geometry_maxRadiusInsideBox(box_width, box_height, (fig_num))
%
% INPUTS:
%      box_width, box_height: the scalar width and height of the box.
%
%      (OPTIONAL INPUTS)
% 
%      None
%
% OUTPUTS:
%
%      max_radius: the radius that fits within the box. Any circles with
%      radii smaller than this will have an arc that does not pass through
%      the sides of the box.
%
% DEPENDENCIES:
%
%      None
%
% EXAMPLES:
%      
% See the script: script_test_fcn_geometry_maxRadiusInsideBox
% for a full test suite.
%
% This function was written on 2024_06_27 by S. Brennan
% Questions or comments? sbrennan@psu.edu 

% Revision history:
% 2024_06_27 - S. Brennan
% -- wrote the code


%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==3 && isequal(varargin{end},-1))
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
        narginchk(2,3);

    end
end


% % Does user want to specify fitting_tolerance?
% fitting_tolerance = 0.1;
% if (2<=nargin)
%     temp = varargin{1};
%     if ~isempty(temp)
%         fitting_tolerance = temp;
%     end
% end
% 
% % Does user want to specify flag_fit_backwards?
% flag_fit_backwards = 0;
% if (3<=nargin)
%     temp = varargin{2};
%     if ~isempty(temp)
%         flag_fit_backwards = temp;
%     end
% end



% Does user want to specify fig_num?
fig_num = []; % Default is to have no figure
flag_do_plots = 0;
if  (0==flag_max_speed) && (3<= nargin)
    temp = varargin{end};
    if ~isempty(temp)        
        fig_num = temp;
        flag_do_plots = 1;
    end
end


%% Solve for the radius
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _       
%  |  \/  |     (_)      
%  | \  / | __ _ _ _ __  
%  | |\/| |/ _` | | '_ \ 
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% METHOD:
% See diagram: https://mathcentral.uregina.ca/QQ/database/QQ.09.07/s/bruce1.html
%
% Let r be the radius, h be the height of the box, and L be the width of
% the box. 
%
% So, the right-triangle in the diagram is formed by a diagonal of length
% r, a base of length 1/2*L, and a height of length (r-h). Using
% Pythagorean theorem:
%
% r^2 = L^2/4 + (r-h)^2
%
% Or:
%
% r^2 = L^2/4 + r^2 - 2*r*h + h^2
%
% The r^2 terms cancel, and rearranging:
%
% 2*r*h = L^2/4 + h^2
%
% or r =(L^2/4 + h^2)/(2*h)
%
% Note: if h>L/2, then the arc can fit inside the box and we need only to
% get the maximum radius, which is just half the box width

L = box_width;
h = box_height;

if h>L/2
    max_radius = L/2;
else
    max_radius = (L^2/4 + h^2)/(2*h);
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

    % Set up figure
    figure(fig_num);
    hold on;
    grid on;
    axis equal;
    xlabel('X [meters]');
    ylabel('Y [meters]');

    % Plot the input square
    square_points = [
        -1*L/2 0
        +1*L/2 0
        +1*L/2 h
        -1*L/2 h
        -1*L/2 0        
        ];
    plot(square_points(:,1),square_points(:,2),'r-','LineWidth',3);

    % Plot the circle
    angles = linspace(0,2*pi,100)';
    circle_center = [0 max_radius];
    circle_points = circle_center + max_radius*[cos(angles) sin(angles)];
    plot(circle_points(:,1),circle_points(:,2),'b-','LineWidth',3);

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