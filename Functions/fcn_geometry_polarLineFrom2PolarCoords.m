function [phi,rho] = fcn_geometry_polarLineFrom2PolarCoords(...
    points,...
    varargin)
% fcn_geometry_polarLineFrom2PolarCoords
% Converts two polar points into the polar form of a line, giving phi and
% rho values for the line.
% See: http://www.nabla.hr/Z_MemoHU-015.htm
%
% FORMAT:
%
% [phi,rho] = fcn_geometry_polarLineFrom2PolarCoords(...
%      points,
%     (fig_num))
%
% INPUTS:
%
%      points: an [2 x 2] vector of [theta1 r1; theta2 r2] which are the
%      ends of the line segment to be fit with a line
%
%      (OPTIONAL INPUTS)
% 
%      fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed.
%
% OUTPUTS:
%
%      [phi,rho] = the phi and rho term in the polar form of a line:
%      cos(theta - phi) = rho/r, where rho is the distance of the line from
%      the origin, and phi is the angle that the perpendicular from the
%      origin to the line makes with the polar axis (e.g. the x-axis)
%
% DEPENDENCIES:
%
%      fcn_DebugTools_checkInputsToFunctions
%      fcn_geometry_plotCircle
%
% EXAMPLES:
%
% See the script: script_test_fcn_geometry_polarLineFrom2PolarCoords
% for a full test suite.
%
% This function was written on 2021_05_22 by S. Brennan
% Questions or comments? sbrennan@psu.edu 

% Revision History:
% 2021-05-27
% -- Created function from fcn_geometry_find_phi_rho_from_two_polar_coords
% 2024_01_17 - Aneesh Batchu
% -- added max speed options



%% Debugging and Input checks
% flag_check_inputs = 1; % Set equal to 1 to check the input arguments
% flag_do_plot = 0;      % Set equal to 1 for plotting
% flag_do_debug = 0;     % Set equal to 1 for debugging

flag_max_speed = 0;
if (nargin==2 && isequal(varargin{end},-1))
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

if 0==flag_max_speed
    if flag_check_inputs
        % Are there the right number of inputs?
        narginchk(1,2);

        % Check the points input, and force it to be a 2x2 matrix
        fcn_DebugTools_checkInputsToFunctions(...
            points, '2column_of_numbers',[2 2]);
    end
end


% % Does user want to show the plots?
% if 2 == nargin
%     fig_num = varargin{end};
%     figure(fig_num);
%     flag_do_plot = 1;
% else
%     if flag_do_debug
%         fig = figure;
%         fig_num = fig.Number;
%         flag_do_plot = 1;
%     end
% end

% Does user want to show the plots?
flag_do_plot = 0;
if (0==flag_max_speed) && (2 == nargin) 
    temp = varargin{1};
    if ~isempty(temp)
        fig_num = temp;
        figure(fig_num);
        flag_do_plot = 1;
    end
else
    if flag_do_debug
        fig = figure; 
        fig_num = fig.Number;
        flag_do_plot = 1;
    end
end

%% Start of main code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _
%  |  \/  |     (_)
%  | \  / | __ _ _ _ __
%  | |\/| |/ _` | | '_ \
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preliminary calculations
theta1 = points(1,1);
r1     = points(1,2);
theta2 = points(2,1);
r2     = points(2,2);

% Fill in the terms
yterms = r1*cos(theta1)-r2*cos(theta2);
xterms = r2*sin(theta2)-r1*sin(theta1);

% Calculate phi and rho
phi = atan2(yterms,xterms);
rho = r1*cos(theta1-phi);


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
if flag_do_plot
    figure(fig_num);
    clf;
    hold on;
    grid on; grid minor;
    
    % Plot the input points
    [Xs,Ys] = pol2cart(points(:,1),points(:,2));
    plot(Xs,Ys,'-o','Markersize',30);
    
    % Fill in a range of angles
    angles = (min(theta1,theta2):0.01:max(theta1,theta2))';
    
    % Use the line formula to calculate the radii for each angle    
    radii_on_segment = rho./cos(angles - phi);
    
    % Convert these polar values into cartesian
    [x_points,y_points] = pol2cart(angles,radii_on_segment);
    
    % Plot the results
    plot(x_points,y_points,'rx');
end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file); 
end

end % Ends function