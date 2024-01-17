function [dist] = ...
    fcn_geometry_euclideanPointsToPointsDistance(...
    points1,...
    points2,...
    varargin)
% fcn_geometry_euclideanPointsToPointsDistance calculates the 
% distance(s) between a vector of points, POINTS1, and another vector of
% points, POINTS2.
%
% FORMAT:
%
% [DIST] = fcn_geometry_euclideanPointsToPointsDistance(POINTS1,POINTS2,(fig_num))
%
% INPUTS:
%
%      POINTS1: an Nx2 or Nx3 series of xy or xyz points 
%      in the form: [x1 y1 z1; x2 y2 z2; ... ; xn yn zn]
%
%      POINTS2: an Nx2 or Nx3 series of xy or xyz points 
%      in the form: [x1 y1 z1; x2 y2 z2; ... ; xn yn zn]
%
%      (OPTIONAL INPUTS)
%
%      fig_num: a figure number to plot results.
%
% OUTPUTS:
%
%      DIST: an N x  1 vector of distances [d1; d2; ... ; dn], where N is
%      the number of point sets
%
% DEPENDENCIES:
%
%      fcn_DebugTools_checkInputsToFunctions
%
% EXAMPLES:
%
%         pt1 = [1 1 5; 5 3 64; 7 2 -2];
%         pt2 = [0 -3 -6; 34 1 17; 18 7 0];
%         dist=fcn_geometry_euclideanPointsToPointsDistance(pt1,pt2);
%
% See the script: script_test_fcn_geometry_euclideanPointsToPointsDistance
% for a full test suite.
%
% This function was written on 2018_11_17 by Seth Tau
% Questions or comments? sat5340@psu.edu 

% Revision History:
% 2021-05-28 - S. Brennan
% -- revised function to prep for geometry class 
% -- rewrote function to use vector sum
% -- added plotting option
% 2021-06-05
% -- fixed comments, added debugging option


%% Debugging and Input checks
flag_check_inputs = 1; % Set equal to 1 to check the input arguments
flag_do_plot = 0;      % Set equal to 1 for plotting
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
    narginchk(2,3);
    
    % Check the points1 input
    fcn_DebugTools_checkInputsToFunctions(...
        points1, '2or3column_of_numbers');
    
    % Use number of rows in points1 to calculate Npoints
    Npoints = length(points1(:,1));
    
    % Check the points2 input, forcing length to match points1
    fcn_DebugTools_checkInputsToFunctions(...
        points2, '2or3column_of_numbers',Npoints);       
end
    

% Does user want to show the plots?
if 3 == nargin
    fig_num = varargin{end};
    figure(fig_num);
    flag_do_plot = 1;
else
    if flag_do_debug
        fig = figure;
        fig_num = fig.Number;
        flag_do_plot = 1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _
%  |  \/  |     (_)
%  | \  / | __ _ _ _ __
%  | |\/| |/ _` | | '_ \
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dist = sum((points1-points2).^2,2).^0.5;

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
    % Set up the figure
    figure(fig_num);
    clf
    hold on;
    grid on; grid minor;
        
    midpoints = (points1+points2)/2;
    for ith_point=1:Npoints
        % 2D plot?
        if length(midpoints(1,:))==2
            % Plot the points
            xdata = [points1(ith_point,1) points2(ith_point,1)];
            ydata = [points1(ith_point,2) points2(ith_point,2)];
            plot(xdata,ydata,'.-','Linewidth',3,'Markersize',20);
            
            % Label the midpoints
            text(midpoints(ith_point,1),midpoints(ith_point,2),sprintf('d - %.1f',dist(ith_point,1)));
        else
            % Plot the points
            xdata = [points1(ith_point,1) points2(ith_point,1)];
            ydata = [points1(ith_point,2) points2(ith_point,2)];
            zdata = [points1(ith_point,3) points2(ith_point,3)];
            plot3(xdata,ydata,zdata,'.-','Linewidth',3,'Markersize',20);

            % Label the midpoints
            text(midpoints(ith_point,1),midpoints(ith_point,2),midpoints(ith_point,3),sprintf('d - %.1f',dist(ith_point,1)));
            
            % Set to 3D view
            view(3);
        end
        
    end
end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file); 
end


end % Ends the function



