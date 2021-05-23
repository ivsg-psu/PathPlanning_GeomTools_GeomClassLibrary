function fcn_geometry_plotCircle(centers,radii,varargin)   

% fcn_geometry_plotCircle -  plots a circle by creating a vector of angles
% spaced 0.01 radians apart, and plotting this as a line around the
% perimeter.
%
% FORMAT:
%
%     fcn_geometry_plotCircle(...
%     centers,...
%     radii,...
%     (fig_num))
%
% INPUTS:
%
%      centers: an [N x 2] vector in [x y] of the points of circle centers
%
%      radii: a [N x 1] vector of the radii of the circles (to avoid
%      calculation time)
%
%      (OPTIONAL INPUTS)
%
%      fig_num: a figure number to plot results.
%
% OUTPUTS:
%
%      (none)
%
% DEPENDENCIES:
%
%      fcn_geometry_checkInputsToFunctions
%
% EXAMPLES:
%
% See the script: script_test_fcn_geometry_plotCircle
% for a full test suite.
%
% This function was written on 2020_05_22 by S. Brennan
% Questions or comments? sbrennan@psu.edu

% Revision History:
% 2021-05-22
% -- new function from fcn_geometry_findAngleUsing3PointsOnCircle


%% Debugging and Input checks
flag_check_inputs = 1; % Set equal to 1 to check the input arguments
flag_do_plot = 0;      % Set equal to 1 for plotting "extra" figures (not used)
flag_do_debug = 0;     % Set equal to 1 for debugging

if flag_do_debug
    st = dbstack; %#ok<*UNRCH>
    fprintf(1,'STARTING function: %s, in file: %s\n',st(1).name,st(1).file);
end

%% check input arguments?
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
    if nargin < 2 || nargin > 3
        error('Incorrect number of input arguments')
    end
    
    % Check the centers input
    fcn_geometry_checkInputsToFunctions(...
        centers, '2column_of_numbers');
    
    % Use number of radii to calculate the number of centers
    num_circles = length(centers(:,1));
    
    % Check the radii input
    fcn_geometry_checkInputsToFunctions(...
        radii, 'column_of_numbers',num_circles);
    
end
    

% Does user want to specify the figure?
flag_new_figure = 0;
if 3 == nargin
    fig_num = varargin{1};
    figure(fig_num);
    flag_do_plot = 1;
else
    if flag_do_debug
        fig = figure;
        fig_num = fig.Number;
        flag_do_plot = 1;
        flag_new_figure = 1;
    else
        fig = gcf;
        fig_num = fig.Number;
    end
end

figure(fig_num);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _
%  |  \/  |     (_)
%  | \  / | __ _ _ _ __
%  | |\/| |/ _` | | '_ \
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if flag_new_figure        
    hold on;
    axis equal
    grid minor;
else
    hold on;
    axis equal
end

angles = 0:0.01:2*pi;
for ith_circle = 1:length(centers(:,1))
    x_circle = centers(ith_circle,1) + radii(ith_circle) * cos(angles);
    y_circle = centers(ith_circle,2) + radii(ith_circle) * sin(angles);
    plot(x_circle,y_circle,'-');
end

%% Plot results?
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
    % Nothing more to do here!
end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file); 
end
