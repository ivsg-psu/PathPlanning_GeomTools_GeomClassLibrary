function [...
    angles] ...
    = ...
    fcn_geometry_findAngleUsing2PointsOnCircle(...
    centers,...
    radii,...
    start_points_on_circle,...
    end_points_on_circle,...
    cross_product_direction,...
    varargin)

% fcn_geometry_findAngleUsing2PointsOnCircle -  This function calculates
% the angle from the start_points location to the end_points, in the
% direction of the vector given by is_clockwise.
%
% FORMAT:
%
% [angles] ...
%     = ...
%     fcn_geometry_findAngleUsing2PointsOnCircle(...
%     centers,...
%     radii,...
%     start_points_on_circle,...
%     end_points_on_circle,...
%     cross_products,...
%     varargin)
%
% INPUTS:
%
%      centers: an [N x 2] vector in [x y] of the points of circle centers
%
%      radii: a [N x 1] vector of the radii of the circles (to avoid
%      calculation time)
%
%      start_points_on_circle: an [N x 2] vector in [x y] of the points
%      where sectors start
%
%      end_points_on_circle: an [N x 2] vector in [x y] of the points
%      where sectors end
%
%      cross_product_direction: an [N x 1] vector denoting the cross product
%      direction to follow from input point to output point
%
%      (OPTIONAL INPUTS)
%
%      fig_num: a figure number to plot results.
%
% OUTPUTS:
%
%      angles: an [N x 1] vector of the angles, in radians, between input
%      points and output points
%
% DEPENDENCIES:
%
%      fcn_geometry_checkInputsToFunctions
%      fcn_geometry_plotCircle
%
% EXAMPLES:
%
% See the script: script_test_fcn_geometry_findAngleUsing2PointsOnCircle
% for a full test suite.
%
% This function was written on 2020_05_22 by S. Brennan
% Questions or comments? sbrennan@psu.edu

% Revision History:
% 2021-05-22
% -- new function from fcn_geometry_findAngleUsing3PointsOnCircle


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
    narginchk(5, 6);
    
    % Check the centers input
    fcn_geometry_checkInputsToFunctions(...
        centers, '2column_of_numbers');
    
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
        cross_product_direction, 'column_of_numbers',num_circles);
       
end
    

% Does user want to show the plots?
if 6 == nargin
    fig_num = varargin{end};
    figure(fig_num);
    flag_do_plot = 1;
    flag_new_figure = 0;
else
    if flag_do_debug
        fig = figure;
        fig_num = fig.Number;
        flag_do_plot = 1;
        flag_new_figure = 1;
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



%% Step 1: calculate unit vectors for incoming and outgoing points
% Makes things easy later
unit_radial_to_inpoints = (start_points_on_circle - centers)./radii;
unit_radial_to_outpoints = (end_points_on_circle - centers)./radii;

%% Step 2: calculate the dot product angle from in and out unit vectors
dot_product = sum(unit_radial_to_inpoints.*unit_radial_to_outpoints,2);

%% Step 3: calculate the cross products from in to out
cross_in_to_out = cross(...
    [unit_radial_to_inpoints, zeros(num_circles,1)],...
    [unit_radial_to_outpoints, zeros(num_circles,1)]);
cross_product = cross_in_to_out(:,3);
cross_product(cross_product==0) = cross_product_direction(cross_product==0);

% Combine them to get the angle
angles = acos(dot_product).*sign(asin(cross_product));

%% Step 4: check if the cross product matches
% If the cross_in_to_out is in opposite direction from
% given cross products, then we need to take the refelex angle.
need_reflex_angles = (cross_product.*cross_product_direction)<0;
angles_positive = angles>0;
angles_negative = ~angles_positive;
negative_reflex_indicies = find(need_reflex_angles.*angles_negative);
positive_reflex_indicies = find(need_reflex_angles.*angles_positive);

angles(negative_reflex_indicies)= 2*pi + angles(negative_reflex_indicies);
angles(positive_reflex_indicies)= angles(positive_reflex_indicies) - 2*pi;

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
        if cross_product_direction(i,1)>0
            text(location(1,1),location(1,2),sprintf(':  %.1f deg, counterclockwise',angles(i,1)*180/pi));
        else
            text(location(1,1),location(1,2),sprintf(':  %.1f deg, clockwise',angles(i,1)*180/pi));
        end
    end
end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file); 
end


end % Ends the function


