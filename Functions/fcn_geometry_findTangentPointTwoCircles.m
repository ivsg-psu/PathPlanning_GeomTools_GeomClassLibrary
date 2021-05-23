function [...
    points_tangent_start, ...
    points_tangent_end] ...
    = ...
    fcn_geometry_findTangentPointTwoCircles(...
    centers_start,...
    centers_end,...
    radii_start,...
    radii_end,...
    cross_products_start,...
    cross_products_end,...
    varargin)
% fcn_geometry_findTangentPointTwoCircles
% finds tangent points from one set of circles to another, returning only
% the one set of points tht matches the given cross products
%
% This function allows vectorization where centers and radii can be a
% vector [x y] and [r] respectively, where x, y, and r are equal-length
% columns. The number of starting centers must match the number of ending
% centers, as the circles are matched to each other when calculating
% outputs. The length of the cross products column vector must also match.
% 
% FORMAT:
%
% [...
%     points_tangent_start, ...
%     points_tangent_end] ...
%     = ...
%     fcn_geometry_findTangentPointTwoCircles(...
%     centers_start,...
%     centers_end,...
%     radii_start,...
%     radii_end,...
%     cross_products_start,...
%     cross_products_end,...
%     (fig_num))
%
% INPUTS:
%
%      centers_start: an [N x 2] vector of X,Y data for each circle center
%      to start the tangent line.
%
%      centers_end: an [N x 2] vector of X,Y data for each circle center
%      to end the tangent line.
%
%      radii_start: a [N x 1] vector of radii for each starting circle
%
%      radii_end: a [N x 1] vector of radii for each starting circle
%
%      cross_products_start: a [N x 1] vector of cross products from the
%      tangent lines (in start-to-end directions) to the circle centers for
%      the start circles. This is used to specify whether there is an inner
%      or outer tangent  and which points to keep. The result is only ONE
%      set of start and end tangent points for each circle: the one that
%      has the matching cross product.
%
%      cross_products_end: a [N x 1] vector of cross products from the
%      tangent lines (in start-to-end directions) to the circle centers for
%      the end circles. This is used to specify whether there is an inner
%      or outer tangent  and which points to keep. The result is only ONE
%      set of start and end tangent points for each circle: the one that
%      has the matching cross product.
%
%      (OPTIONAL INPUTS)
%  
%      fig_num: a figure number to plot results.
%
% OUTPUTS:
%
%      points_tangent_start: an [N x 2] vector of X,Y data for each circle.
%
%      points_tangent_end: an [N x 2] vector of X,Y data for each circle.
%
%
% DEPENDENCIES:
%
%      fcn_geometry_checkInputsToFunctions
%      fcn_geometry_findTangentPointsFromPointToCircle
%      fcn_geometry_plotCircle
%
% EXAMPLES:
%
% See the script: script_test_fcn_geometry_findTangentPointTwoCircles
% for a full test suite.
%
% This function was written on 2021_04_25 by S. Brennan
% by modifying: fcn_geometry_findTangentPointsTwoCircles
% Questions or comments? sbrennan@psu.edu 
%
% See: http://www.ambrsoft.com/TrigoCalc/Circles2/Circles2Tangent_.htm for
% the mathematics being implemented here.

% Revision History:
% 2021-04-25
% -- First write of the code
% 2021-05-22
% -- Added plotting from: fcn_geometry_plotCircle

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
% The input arguments are:
% (centers_start,centers_end,radii_start,radii_end,flag_inside_or_outside,voting_points_start,voting_points_end,fignum);

Ncircles = length(centers_start(:,1));  % The number of start and end circles

if flag_check_inputs
    if nargin < 6 || nargin > 7
        error('Incorrect number of input arguments.')
    end
    
    % Check the centers_start input
    fcn_geometry_checkInputsToFunctions(...
        centers_start, '2column_of_numbers',Ncircles);
    
    % Check the centers_end input
    fcn_geometry_checkInputsToFunctions(...
        centers_end, '2column_of_numbers',Ncircles);

    % Check the radii_start input
    fcn_geometry_checkInputsToFunctions(...
        radii_start, 'column_of_numbers',Ncircles);

    % Check the radii_end input
    fcn_geometry_checkInputsToFunctions(...
        radii_end, 'column_of_numbers',Ncircles);
    
    % Check the cross_products_start input
    fcn_geometry_checkInputsToFunctions(...
        cross_products_start, 'column_of_numbers',Ncircles);
        
    % Check the cross_products_end input
    fcn_geometry_checkInputsToFunctions(...
        cross_products_end, 'column_of_numbers',Ncircles);
    
    
end

% Does user want to show the plots?
if 7 == nargin
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
%% Initialize vectors
points_tangent_start = nan(Ncircles,2);
points_tangent_end   = nan(Ncircles,2);

% Find situations that are inside or outside
inside_or_outside = cross_products_start.*cross_products_end;

% Check if any zeros
bad_indices = find(inside_or_outside==0, 1);
if ~isempty(bad_indices)
    % This only happens if one of the cross products is zero. And this only
    % happens if a path is given that does not have a bend, or bends
    % perfectly backwards. In either case, it is impossible to determine
    % whether to use an inside or outside tangent circle and so an error is thrown.
    error('Trying to calculate tangents to circles where the cross product is zero. This occurs when one portion of the path has no bend at an apex, or is bending backwards. Impossible to determine how to create internal or external tangents without cross products.');   
end

inside_indices  = find(inside_or_outside<0);
outside_indices = find(inside_or_outside>0);

%% Calculate results for inside intersection points

if ~isempty(inside_indices)
    % Check that circles are far enough apart to have inside tangents
    r_sum_squared = ...
        (radii_start(inside_indices,:) +...
        radii_end(inside_indices,:)).^2;
    D_squared = sum(...
        (centers_start(inside_indices,:) - ...
        centers_end(inside_indices,:)).^2,2);
    if any(D_squared<r_sum_squared)
        % Plot the circles
        fcn_plotCircles(centers_start,centers_end,radii_start,radii_end,fig_num);
        error('The internal tangent points between two circles were sought, but the circles overlap each another. Not possible to continue!');
    end
    
    % Calculate inner tangent points
    points_inner = ...
        (centers_end.*radii_start + ...
        centers_start.*radii_end)./ ...
        (radii_start + radii_end);
    
    % Find intersection points using a point and circle formula. Since the
    % inner tanget point is between start and end, and cross product in the
    % function below is calculated from the point to the circle, the cross
    % product will have the opposite sign than normal.
    points_tangent_start(inside_indices,:) = ...
        fcn_geometry_findTangentPointFromPointToCircle(...
        centers_start(inside_indices,:),...
        radii_start(inside_indices,:),...
        points_inner(inside_indices,:),...
        -cross_products_start(inside_indices,:));
    
    % Do the calculation for the end circle. In this case, the cross
    % product is in the correct direction
    points_tangent_end(inside_indices,:)   = ...      
        fcn_geometry_findTangentPointFromPointToCircle(...
        centers_end(inside_indices,:),...
        radii_end(inside_indices,:),...
        points_inner(inside_indices,:),...
        cross_products_end(inside_indices,:));   
      
end % Ends check if inside tangent points should be calculated

%% Calculate outer tangent intersection points. 
if ~isempty(outside_indices)
    

    % Check to see if these are equi-radii situations. If the test circles
    % have the same radius, then the methods below won't work. So here, we
    % perform special calculations for circles of the same radius, only for
    % data where the radii are equal (noted by the variable equal_indices).
    
    % First, check the radii difference from start and end circle pairs
    radii_difference = radii_start - radii_end;
   
    % flag if radii difference is less than tolerance.
    tol = 0.00001;
    equal_indices = find((radii_difference.^2)<tol^2);
    
    % Keep only the equal indices that are also outside_indices
    equal_and_outside_indices = intersect(outside_indices,equal_indices);
    unequal_and_outside_indices = setdiff(outside_indices,equal_indices);
   
    % If any flags exist, do calculations (special) for same circle
    % situations
    if any(equal_and_outside_indices)
        
        % Calculate the angle of the line connecting the center of the
        % circles
        diff_centers = ...
            centers_end(equal_and_outside_indices,:) - ...
            centers_start(equal_and_outside_indices,:);
        angles_centers = atan2(diff_centers(:,2),diff_centers(:,1));
        
        % Create 2 new angles, each +/- 90 degrees from that circle-circle
        % angle
        angles_tangents1 = angles_centers+pi/2; % This gives the negative cross product
        angles_tangents2 = angles_centers-pi/2; % This gives the positive cross product
        
        % Now calculate where the tangent points would be, using the 2 angles
        % above, using a polar projection        
        points_tangent_start(equal_and_outside_indices,:) =...
            centers_start(equal_and_outside_indices,:) + ...
            radii_start(equal_and_outside_indices,:).*(...
            [cos(angles_tangents1) sin(angles_tangents1)].*...
            (cross_products_start(equal_and_outside_indices,:)<0) + ...
            [cos(angles_tangents2) sin(angles_tangents2)].*...
            (cross_products_start(equal_and_outside_indices,:)>=0));
        
        
        
        points_tangent_end(equal_and_outside_indices,:) =...
            centers_end(equal_and_outside_indices,:) + ...
            radii_end(equal_and_outside_indices,:).*(...
            [cos(angles_tangents1) sin(angles_tangents1)].*...
            (cross_products_end(equal_and_outside_indices,:)<0) + ...
            [cos(angles_tangents2) sin(angles_tangents2)].*...
            (cross_products_end(equal_and_outside_indices,:)>=0));
        
           
    end % Ends outside tangents for equal radii
    
    % Calculate outer tangent points for non-equal-radii case (more typical)
    if ~isempty(unequal_and_outside_indices)
        % Check to see if one circle is inside the other
        r_diff_squared = radii_difference(unequal_and_outside_indices,:).^2;
        D_squared = sum(...
            (centers_start(unequal_and_outside_indices,:) - ...
            centers_end(unequal_and_outside_indices,:)).^2,2);
        if any(D_squared<r_diff_squared)
            % Plot the circles
            fcn_plotCircles(centers_start,centers_end,radii_start,radii_end,fig_num);
            error('The external tangent points between two circles were sought, but one circle is completely inside another. Not possible to continue!');
        end
               
        % If no circles are within another, now can find the outer points
        points_outer = ...
            (centers_end(unequal_and_outside_indices,:).*...
            radii_start(unequal_and_outside_indices,1) - ...
            centers_start(unequal_and_outside_indices,:).*...
            radii_end(unequal_and_outside_indices,1))./...
            radii_difference(unequal_and_outside_indices,1);
        
        % The cross products must switch signs depending on whether the
        % start circle is greater in radius than the end circle
        cross_product_start_modified = ...
            cross_products_start(unequal_and_outside_indices,:).* ...
            ((radii_difference(unequal_and_outside_indices,:)<0) - ...
            (radii_difference(unequal_and_outside_indices,:)>0));
        cross_product_end_modified = ...
            cross_products_start(unequal_and_outside_indices,:).* ...
            ((radii_difference(unequal_and_outside_indices,:)<0) - ...
            (radii_difference(unequal_and_outside_indices,:)>0));
        
        % Find intersection points using a point and circle tangent, using
        % the fixed cross-products to down-select
         points_tangent_start(unequal_and_outside_indices,:) = ...
             fcn_geometry_findTangentPointFromPointToCircle(...
             centers_start(unequal_and_outside_indices,:),...
             radii_start(unequal_and_outside_indices,:),...
             points_outer,...
             cross_product_start_modified);        
         points_tangent_end(unequal_and_outside_indices,:) = ...
             fcn_geometry_findTangentPointFromPointToCircle(...
             centers_end(unequal_and_outside_indices,:),...
             radii_end(unequal_and_outside_indices,:),...
             points_outer,...
             cross_product_end_modified);
    end % Ends outside tangents, unequal radii
    
end % Ends check to see if outside tangents should be calculated


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
    hold on;
    axis equal;
    grid on; grid minor;
    
    % Plot the circles
    fcn_plotCircles(centers_start,centers_end,radii_start,radii_end,fig_num);
    
    % Plot the tangent points    
    plot(points_tangent_start(:,1),points_tangent_start(:,2),'g+');
    plot(points_tangent_end(:,1),points_tangent_end(:,2),'g+');
    
    % Plot the tangent lines
    for  i=1:length(points_tangent_start(:,1))
        plot([points_tangent_start(i,1) points_tangent_end(i,1)],...
            [points_tangent_start(i,2) points_tangent_end(i,2)],...
            'r-');
    end
    
end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file); 
end

end % Ends function

%% For plotting
function fcn_plotCircles(centers_start,centers_end,radii_start,radii_end,fig_num)
figure(fig_num);
hold on;
axis equal;
grid on; grid minor;

% Plot the circle centers
plot(centers_start(:,1),centers_start(:,2),'+');
plot(centers_end(:,1),centers_end(:,2),'+');

% plot the circles, and label then with S for start, E for end
fcn_geometry_plotCircle(centers_start,radii_start);
text(centers_start(:,1),centers_start(:,2)+radii_start*0.5,'S');
fcn_geometry_plotCircle(centers_end,radii_end);
text(centers_end(:,1),centers_end(:,2)-radii_start*0.5,'E');
end



