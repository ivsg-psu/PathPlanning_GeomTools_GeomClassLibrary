function fittedSpiral = fcn_geometry_fitSpiralBetween2points(endPoint, startPoint, varargin)
% fcn_geometry_fitArc
%
% This function takes the end point of a/an straight line/arc and the start
% point of an/a arc/straight line as the input and outputs the fitted
% points of the Archimedean spiral. In short, this function fits an
% Archimedean spiral to the end point of a/an straight line/arc and start
% point of an arc/straight line
%
% METHOD
%
% STEP 1: Determine the angles of the end and start points
%
% STEP 2: Compute a (Initial distance) and b (Rate of growth) parameter by
% calculating the distances to the end and start points from the origin.
%
% STEP 3: Using a and b, fit an Archimedean spiral (a + b*theta) b/w the
% end and start points
%
%
% FORMAT:
%
%       fittedSpiral = fcn_geometry_fitSpiralBetween2points(endPoint, startPoint, varargin)
%
% INPUTS: 
%       endPoint: The end point of a/an straight line/arc
%
%       startPoint: The start point of an/a arc/straight line
%
%
% (OPTIONAL INPUTS)
%
%       fig_num: The figure is plotted if fig_num is entered as the input. 
%
% OUTPUTS:
% 
%       fittedSpiral: The fitted points of an Archimedean Spiral
%
%
% DEPENDENCIES: 
% 
%       (None)
% 
% EXAMPLES:
%
%     See the script: script_test_fcn_Transform_ENUToSensorCoord
%     for a full test suite.

% REVISION HISTORY
% 
% 2023_10_19: Aneesh Batchu
% -- wrote this code originally


flag_do_debug = 0;  % Flag to show the results for debugging
flag_do_plots = 0;  % % Flag to plot the final results
flag_check_inputs = 1; % Flag to perform input checking
fileID = 1; % The default file ID destination for fprintf messages

if flag_do_debug
    st = dbstack; %#ok<*UNRCH>
    fprintf(fileID,'STARTING function: %s, in file: %s\n',st(1).name,st(1).file);
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

if flag_check_inputs
    % Are there the right number of inputs?
    narginchk(2,3);

end

% Does user want to show the plots?
% fig_num = [];
if 3 == nargin
    temp = varargin{end};
    if ~isempty(temp)
        fig_num = temp;
        figure(fig_num);
        flag_do_plots = 1; 
    end
else
    if flag_do_debug
        flag_do_plots = 1;
    end
end

%% Main code starts here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _
%  |  \/  |     (_)
%  | \  / | __ _ _ _ __
%  | |\/| |/ _` | | '_ \
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate the angles of the points
thetaEndPoint = atan2(endPoint(1,2), endPoint(1,1));
thetaStartPoint = atan2(startPoint(1,2), startPoint(1,1));

% Calculate the values of a and b
r_endPoint = sqrt(sum(endPoint.^2));
r_startPoint = sqrt(sum(startPoint.^2));
b = (r_startPoint - r_endPoint) / (thetaStartPoint - thetaEndPoint);
a = r_endPoint - b * thetaEndPoint;

% Generate points along the fitted Archimedean spiral
theta = linspace(thetaEndPoint, thetaStartPoint, 50);
r = a + b * theta;
x = r .* cos(theta);
y = r .* sin(theta);

fittedSpiral = [x' y'];


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

    hold on
    grid on
    box on
    axis equal

    title('Fitted Archimedean Spiral between Two Points');

    % Plot the fitted Archimedean spiral and the two points
    plot([endPoint(1,1), startPoint(1,1)], [endPoint(1,2), startPoint(1,2)], 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g');

    plot(fittedSpiral(:,1), fittedSpiral(:,2), 'r', 'LineWidth', 2);

end

if flag_do_debug
    fprintf(fileID,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file);
end

end

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