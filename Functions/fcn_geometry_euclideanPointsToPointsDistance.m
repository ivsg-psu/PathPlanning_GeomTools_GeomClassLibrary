function [dist] = fcn_geometry_euclideanPointsToPointsDistance(pt1,pt2, varargin)
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
%      fcn_geometry_checkInputsToFunctions
%
% EXAMPLES:
%
%         pt1 = [1 1 5; 5 3 64; 7 2 -2];
%         pt2 = [0 -3 -6; 34 1 17; 18 7 0];
%         dist=fcn_geometry_euclideanPointsToPointsDistance(pt1,pt2);
%
% See the script: script_test_fcn_geometry_findAngleUsing2PointsOnCircle
% for a full test suite.
%
% This function was written on 2018_11_17 by Seth Tau
% Questions or comments? sat5340@psu.edu 

% Revision History:
% 2021-05-28 - S. Brennan
% -- revised function to prep for geometry class 
% -- rewrote function to use vector sum
% -- added plotting option


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
    if nargin < 2 || nargin > 3
        error('Incorrect number of input arguments')
    end
    
    % Check the pt1 input
    fcn_geometry_checkInputsToFunctions(...
        pt1, '2or3column_of_numbers');
    
    % Use number of radii to calculate the number of centers
    num_circles = length(centers(:,1));
    
    % Check the radii input
    fcn_geometry_checkInputsToFunctions(...
        radii, 'column_of_numbers',num_circles);
    
    % Check the start_points_on_circle input
    fcn_geometry_checkInputsToFunctions(...
        start_points_on_circle, '2column_of_numbers',num_circles);
    
    % Check the end_points_on_circle input
    fcn_geometry_checkInputsToFunctions(...
        end_points_on_circle, '2column_of_numbers',num_circles);
        
    % Check the cross_products input
    fcn_geometry_checkInputsToFunctions(...
        cross_products, 'column_of_numbers',num_circles);
       
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

if col1 == 2
    dist = sqrt((pt1(:,1)-pt2(:,1)).^2 + (pt1(:,2)-pt2(:,2)).^2);
else
    dist = sqrt((pt1(:,1)-pt2(:,1)).^2 + (pt1(:,2)-pt2(:,2)).^2 + (pt1(:,3)-pt2(:,3)).^2);
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
    % Set up the figure
    if flag_new_figure
        figure(fig_num);
    else
        clf;
    end
    hold on;
    axis equal;
    grid on; grid minor;
    
    
    % Plot the circles
    fcn_geometry_plotCircle(centers,radii);
    axis equal;
    
    plot(centers(:,1),centers(:,2),'kx');
    
    % Plot the start and end points
    plot(start_points_on_circle(:,1),start_points_on_circle(:,2),'go');
    text(start_points_on_circle(:,1),start_points_on_circle(:,2),'Start');
    
    plot(end_points_on_circle(:,1),end_points_on_circle(:,2),'rx');
    text(end_points_on_circle(:,1),end_points_on_circle(:,2),'End');
    
    for i=1:num_circles
        % Plot unit vectors
        plot([centers(i,1); centers(i,1)+radii(i).*unit_radial_to_inpoints(i,1)],...
            [centers(i,2); centers(i,2)+radii(i).*unit_radial_to_inpoints(i,2)],'g');
        
        plot([centers(i,1); centers(i,1)+radii(i).*unit_radial_to_outpoints(i,1)],...
            [centers(i,2); centers(i,2)+radii(i).*unit_radial_to_outpoints(i,2)],'r');
        
        % Plot the angle value
        location = centers(i,:);
        if cross_products(i,1)>0
            text(location(1,1),location(1,2),sprintf(':  %.1f deg counterclockwise',angles(i,1)*180/pi));
        else
            text(location(1,1),location(1,2),sprintf(':  %.1f deg clockwise',angles(i,1)*180/pi));
        end
    end
end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file); 
end


end % Ends the function



