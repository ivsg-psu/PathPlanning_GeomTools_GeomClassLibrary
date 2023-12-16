function single_points_tangent = ...
    fcn_geometry_findTangentPointFromPointToCircle(...
    centers,...
    radii,...
    points,...
    cross_product_goal,...
    varargin)
% fcn_geometry_findTangentPointFromPointToCircle
% finds the ONE tangent point on a circle given a point and a cross product.
%
% This calculates the tangent points on circles defined by centers and
% radii, where tangents pass through location given by points, and keeps
% only the point that has the same sign as the given cross product. This
% function allows vectorization where centers and points can be a vector [x
% y] where x and y are columns. The radii must be a column vector of same
% length as centers. Points are also in [x y] format, of same length as
% centers. 
% 
% FORMAT:
%
% points_tangent = ...
%     fcn_geometry_findTangentPointFromPointToCircle(...
%     centers,...
%     radii,...
%     points,...
%     cross_product_goal,...
%     (fig_num))
%
% INPUTS:
%
%      centers: an [N x 2] vector of X,Y data for each circle center
%
%      radii: a [N x 1] vector of radii for each circle center
%
%      points: an [N x 2] vector of X,Y data for each point
%
%      cross_product_goal: a [N x 1] vector of cross products for each circle center
%
%      (OPTIONAL INPUTS)
%
%      fig_num: a figure number to plot results.
%
% OUTPUTS:
%
%      points_tangent: an [N x 2] vector of X,Y data for each circle.
%
% DEPENDENCIES:
%
%      fcn_geometry_checkInputsToFunctions
%      fcn_geometry_findTangentPointsFromPointToCircle
%
% EXAMPLES:
%
% See the script: script_test_fcn_geometry_findTangentPointFromPointToCircle
% for a full test suite.
%
% This function was written on 2020_03_28 by S. Brennan
% Questions or comments? sbrennan@psu.edu 
%
% See: http://www.ambrsoft.com/TrigoCalc/Circles2/Circles2Tangent_.htm for
% the mathematics being implemented here.

% Revision History:
% 2021-04-23
% -- Revised the comments area, prepped function for geometry class
% 2021-04-24
% -- Used the results from fcn_geometry_findTangentPointsFromPointToCircle
% to simplify the code (and force error-checking into that part).



%% Debugging and Input checks
flag_check_inputs = 1; % Set equal to 1 to check the input arguments
flag_do_plot = 0;      % Set equal to 1 for plotting
flag_do_debug = 0;     % Set equal to 1 for debugging

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
if flag_check_inputs    
    % Are there the right number of inputs?
    narginchk(4,5);
    
    % Check the centers input
    fcn_geometry_checkInputsToFunctions(...
        centers, '2column_of_numbers');
    
    Ncircles = length(centers(:,1));
    
    % Check the radii input
    fcn_geometry_checkInputsToFunctions(...
        radii, 'column_of_numbers',Ncircles);
    
    % Check the points input    
    fcn_geometry_checkInputsToFunctions(...
        points, '2column_of_numbers',Ncircles);

    % Check the cross_product_goal input
    fcn_geometry_checkInputsToFunctions(...
        cross_product_goal, 'column_of_numbers',Ncircles);
end


% Does user want to show the plots?
if 5 == nargin
    fig_num = varargin{1};
    figure(fig_num);
    flag_do_plot = 1;
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

% Calculate both the tangent points
points_tangent = ...
    fcn_geometry_findTangentPointsFromPointToCircle(...
    centers,...
    radii,...
    points);

%% Separate out the points
points_tanget1 = points_tangent(1:Ncircles,:);
points_tanget2 = points_tangent(Ncircles+1:end,:);


%% Calculate the vectors and cross products
% Calculate the vectors from the tangent points to the center of the
% circles
vectors_tangent_to_circleCenters_1 = centers - points_tanget1;
vectors_points_to_tangents_1 = points_tangent(1:Ncircles,:) - points;

% Do the cross products
cross_product_tangents_to_circleCenters = cross(...
    [vectors_points_to_tangents_1 zeros(Ncircles,1)],...
    [vectors_tangent_to_circleCenters_1 zeros(Ncircles,1)]);

% Check the signs
signs = (cross_product_tangents_to_circleCenters(:,3).*cross_product_goal)>0;

% Assign the results
single_points_tangent = points_tanget1.*signs + points_tanget2.*(~signs);

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
    axis equal;
    grid on; grid minor;
    plot_circle(centers,radii);
    plot(centers(:,1),centers(:,2),'kx');
    plot(points(:,1),points(:,2),'go');
    plot(points_tangent(:,1),points_tangent(:,2),'gx');
    
    for i=1:Ncircles
         plot([points(i,1) single_points_tangent(i,1)],[points(i,2) single_points_tangent(i,2)],'r-');
         %          plot([points(i,1) points_tangent(i,1)],[points(i,2) points_tangent(i,2)],'r-');
         %          plot([points(i,1) points_tangent(i+Ncircles,1)],[points(i,2) points_tangent(i+Ncircles,2)],'r-');
    end
end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file); 
end

end % Ends function

function plot_circle(centers,radii)    % plot all the circle fits
    angles = 0:0.01:2*pi;
    for i_fit = 1:length(centers(:,1))
        x_circle = centers(i_fit,1) + radii(i_fit) * cos(angles);
        y_circle = centers(i_fit,2) + radii(i_fit) * sin(angles);
        plot(x_circle,y_circle,'b-');
    end
end

