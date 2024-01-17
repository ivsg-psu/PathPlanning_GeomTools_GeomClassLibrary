function [fittedLine, fittedParameters] = fcn_geometry_fitStraightLine(noisyData, varargin)
% fcn_geometry_fitStraightLine
%
% This function takes noisy data as the input and outputs the fitted points
% of the straight line. In short, this function fits a straight line to the
% noisy points. 
%
% METHOD
%
% STEP 1: The polyfit function is used to fit the data and determine the
% parameters to fit a line.
%
% STEP 2: Linspace function is used to generate x-values for the line
%
% STEP 3: polyval is used to compute corresponding y-values for the line
%
% FORMAT:
%
%       [fittedLine, fittedParameters] = fcn_geometry_fitStraightLine(noisyData, varargin)
%
% INPUTS: 
%       noisyData: noisy data to fit the line
%
% (OPTIONAL INPUTS)
%
% fig_num: The figure is plotted if fig_num is entered as the input. 
%
% OUTPUTS:
% 
%       fittedLine: The fitted points of a straight line
%
%       fittedParameters: The parameters that are used to fit a straight
%       line to the noisy points
%       
% DEPENDENCIES: 
% 
%       (None)
% 
% EXAMPLES:
%
%     See the script: script_test_fcn_geometry_fitStraightLine
%     for a full test suite.

% REVISION HISTORY
% 
% 2023_10_16: Aneesh Batchu
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
    narginchk(1,2);

end

% Does user want to show the plots?
% fig_num = [];
if 2 == nargin
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

%Fit a line (y = mx + b) to the data
fittedParameters = polyfit(noisyData(:,1), noisyData(:,2), 1);

% Generate x-values for the line
xFit = linspace(min(noisyData(:,1)), max(noisyData(:,1)), size(noisyData,1));

% Compute corresponding y-values for the line
yFit = polyval(fittedParameters, xFit); 

fittedLine = [xFit', yFit']; 


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
    title('Random Data Fitted to a Straight Line');
    
    % Plot the noisy points
    scatter(noisyData(:,1), noisyData(:,2), 'blue', 'filled');

    % Plot the fitted straight line
    plot(fittedLine(:,1), fittedLine(:,2), 'r', 'LineWidth', 2);
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