function [fittedArc,fittedParameters] = fcn_geometry_fitArc(noisyData, initialGuess, varargin)
% fcn_geometry_fitArc
%
% This function takes noisy data as the input and outputs the fitted points
% of the arc. In short, this function fits an arc to the noisy points. 
%
% METHOD
%
% STEP 1: Determine the angles between x and y oordinates of the noisy
% points using atan2
%
% STEP 2: Create a model function (fcn_arcFit) function to fit an arc
%
% STEP 3: Determine the fitted parameters by calling the lsqcurvefit
% function and by passing model function, initial guess, angles and noisy
% points
%
% STEP 4: Divide the fitted startAngle (parameter(1,4)) and fitted endAngle
% (parameter(1,5)) into size(noisyPoints,1) elements.
%
%
% FORMAT:
%
%       [fittedArc,fittedParameters] = fcn_geometry_fitArc(noisyData, initialGuess, varargin)
%
% INPUTS: 
%       noisyData: noisy data to fit the line
%
%       initialGuess: [initial_guess_xCoor_of_center,
%                      initial_guess_yCoor_of_center,
%                      initial_guess_radius,
%                      initial_guess_start_angle,
%                      initial_guess_end_angle]s
%
% (OPTIONAL INPUTS)
%
%       fig_num: The figure is plotted if fig_num is entered as the input. 
%
% OUTPUTS:
% 
%       fittedArc: The fitted points of an arc
%
%       fittedParameters: The parameters that are used to fit an arc to
%       the noisy points
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
% 2023_10_17: Aneesh Batchu
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

% Generate random angles between startAngle and endAngle
% angles = startAngle + (endAngle - startAngle)  * rand(size(noisyData,1), 1);
angles = atan2(noisyData(:,2) - initialGuess(1,2), noisyData(:,1) - initialGuess(1,1));

% Call lsqcurvefit function to fit an arc to the noisy points
fittedParameters = lsqcurvefit(@fcn_arcFit, initialGuess, angles, noisyData);

% The fitted angles to fit an arc
fittedAngles = linspace(fittedParameters(1,4), fittedParameters(1,5), size(noisyData,1));

% The fitted points of an arc
fittedArc = fcn_arcFit(fittedParameters, fittedAngles');

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

    title('Random Data Fitted to an Arc');
    legend('Noisy Points', 'Fitted Arc', 'Location', 'best');
    
    % Plot the noisy points
    scatter(noisyData(:,1), noisyData(:,2), 'b', 'filled');

    % Plot the fitted arc
    plot(fittedArc(:,1), fittedArc(:,2), 'r', 'LineWidth', 2);

    % Mark the center of the fitted arc
    plot(fittedParameters(1,1), fittedParameters(1,2), 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');

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

% The function of an arc.
function arcFit = fcn_arcFit(parameters, angles)
% Arc Model parameters
% parameters(1,1): X-coordinate of the center
% parameters(1,2): Y-coordinate of the center
% parameters(1,3): Radius of the arc
% parameters(1,4): Start Angle in rad
% parameters(1,5): End Angle in rad

theta = linspace(parameters(1,4), parameters(1,5), size(angles,1));
theta = theta';
xFit_points = parameters(1,1) + parameters(1,3) * cos(theta);
yFit_points = parameters(1,2) + parameters(1,3) * sin(theta);
arcFit = [xFit_points, yFit_points];
end