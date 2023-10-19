function [fittedSpiral,fittedParameters] = fcn_geometry_fitSpiral(noisyData, initialGuess, varargin)
% fcn_geometry_fitSpiral
%
% This function takes noisy data as the input and outputs the fitted points
% of the spiral. In short, this function fits an Archimedian spiral to the
% noisy points.
%
% METHOD
%
% STEP 1: Determine the angles between x and y oordinates of the noisy
% points using atan2
%
% STEP 2: Create a model function (fcn_spiralFit) function to fit a spiral
%
% STEP 3: Determine the fitted parameters by calling the lsqcurvefit
% function and by passing model function, initial guess, angles and noisy
% points
%
% STEP 4: Divide the fitted startAngle (thetaFirstPoint) and fitted endAngle
% (thetaLastPoint) into size(noisyPoints,1) elements.
%
%
% FORMAT:
%
%       [fittedSpiral,fittedParameters] = fcn_geometry_fitSpiral(noisyData, initialGuess, varargin)
%
% INPUTS: 
%       noisyData: noisy data to fit the line
%
%       initialGuess: [initial_guess_distance_a,
%                      initial_guess_rate_of_growth_b]
%
% (OPTIONAL INPUTS)
%
%       fig_num: The figure is plotted if fig_num is entered as the input. 
%
% OUTPUTS:
% 
%       fittedSpiral: The fitted points of an Archimedean spiral
%
%       fittedParameters: The parameters that are used to fit an
%       Archimedean spiral to the noisy points
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


% The angles of the first and the last points
thetaFirstPoint = atan2(noisyData(1,2), noisyData(1,1));
thetaLastPoint = atan2(noisyData(end,2), noisyData(end,1));

% The angles of noisy data
noisyAngles = atan2(noisyData(:,2), noisyData(:,1));

% Call lsqcurvefit function to fit an Archimedean spiral (a + b*theta) to
% the noisy points
fittedParameters = lsqcurvefit(@fcn_spiralFit, initialGuess, noisyAngles, noisyData);

% The fitted angles to fit an Archimedean spiral (a + b*theta)
fittedAngles = linspace(thetaFirstPoint, thetaLastPoint, size(noisyData,1));

% The fitted points of an Archimedean spiral (a + b*theta)
fittedSpiral = fcn_spiralFit(fittedParameters, fittedAngles');

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

    title('Random Data Fitted to an Archimedean spiral');
    
    % Plot the noisy points
    scatter(noisyData(:,1), noisyData(:,2), 'b', 'filled');

    % Plot the fitted Archimedean spiral (a + b*theta)
    plot(fittedSpiral(:,1), fittedSpiral(:,2), 'r', 'LineWidth', 2);

    % Mark the first and the last points of the fitted Archimedean spiral
    plot([noisyData(1,1), noisyData(end,1)] , [noisyData(1,2), noisyData(end,2)], 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g');

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

% The function of an Archimedean spiral (a + b*theta). 
function spiralFit = fcn_spiralFit(parameters, theta)
    % Archimedean spiral model function
    % params(1): Initial distance
    % params(2): Radius scaling factor/ Rate of growth
    r = parameters(1) + parameters(2) * theta;
    xFit = r .* cos(theta);
    yFit = r .* sin(theta);
    spiralFit = [xFit, yFit];
end
