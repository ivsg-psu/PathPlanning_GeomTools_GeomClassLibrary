function [...
    points_tangent_start, ...
    points_tangent_end,...
    inconsistent_votes] ...
    = ...
    fcn_geometry_findTangentPointsTwoCircles(...
    centers_start,...
    centers_end,...
    radii_start,...
    radii_end,...
    varargin)
% fcn_geometry_findTangentPointsTwoCircles
% finds tangent points from one set of circles to another
%
% This calculates all the tangent points on circles defined by centers and
% radii, where tangents can be either inner tangents or outer tangents. An
% optional flag can force one or the other to be used.
%
% This function allows vectorization where centers and radii can be a
% vector [x y] and [r] respectively, where x, y, and r are equal-length
% columns. The number of starting centers must match the number of ending
% centers, as the circles are matched to each other when calculating
% outputs.
% 
% FORMAT:
%
% [...
%     points_tangent_start, ...
%     points_tangent_end] ...
%     = ...
%     fcn_geometry_findTangentPointsTwoCircles(...
%     centers_start,...
%     centers_end,...
%     radii_start,...
%     radii_end,...
%     (flag_inside_or_out),...
%     (voting_points_start,voting_points_end),...
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
%      (OPTIONAL INPUTS)
%
%      flag_inside_or_outside: a scalar flag indicating whether to
%      calculate all the tangent points, or some subset. The flag can be
%      one of the following:
%            flag = 0 (default) to calculate both inside and outside
%            tangents
% 
%            flag = 1 to calculate outside only, 
%
%            flag = -1 to calculate inside tangents only
%
%      voting_points_start,voting_points_end: each is a [N x 2] vector of
%      X,Y points which allows user to enter voting points, to keep only
%      tangents whose start and end are closest to the voting points.
%  
%      fig_num: a figure number to plot results.
%
% OUTPUTS:
%
%      points_tangent_start: an [N x 2] vector of X,Y data for each circle.
%
%      points_tangent_end: an [N x 2] vector of X,Y data for each circle.
%
%      inconsistent_votes: shows consistency check for the votes as an
%      array the same length of votes. Inconsistent votes will be flagged
%      as 1.
%
%
% DEPENDENCIES:
%
%      fcn_DebugTools_checkInputsToFunctions
%      fcn_geometry_findTangentPointsFromPointToCircle
%      fcn_geometry_plotCircle
%
% EXAMPLES:
%
% See the script: script_test_fcn_geometry_findTangentPointsTwoCircles
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
% 2021-05-22
% -- Added plotting from: fcn_geometry_plotCircle
% -- Fixed typo on radii difference

%% Debugging and Input checks
flag_check_inputs = 1; % Set equal to 1 to check the input arguments
flag_do_plot = 0;      % Set equal to 1 for plotting
flag_do_debug = 0;     % Set equal to 1 for debugging
flag_do_voting = 0;    % Flag to do the voting. Overwritten below if that is set
flag_inside_or_outside = 0;  % Flag to do inside or outside calculations. Overwritten below if set by user
voting_points_start = [];
voting_points_end = [];

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
    if nargin < 4 || nargin > 8
        error('Incorrect number of input arguments.')
    end
    
    % Check the centers_start input
    fcn_DebugTools_checkInputsToFunctions(...
        centers_start, '2column_of_numbers',Ncircles);
    
    % Check the centers_end input
    fcn_DebugTools_checkInputsToFunctions(...
        centers_end, '2column_of_numbers',Ncircles);

    % Check the radii_start input
    fcn_DebugTools_checkInputsToFunctions(...
        radii_start, 'column_of_numbers',Ncircles);

    % Check the radii_end input
    fcn_DebugTools_checkInputsToFunctions(...
        radii_end, 'column_of_numbers',Ncircles);
end


if 5 <= nargin  % Flag for inside or outside is set by user
    flag_inside_or_outside = varargin{1};
end

if 6 <= nargin  % User is providing voting points
    if 7>nargin  % user did not give voting points end?
        error('If starting votes vector is specified, end points vector must also be specified');
    end
    voting_points_start = varargin{2};
    voting_points_end = varargin{3};
    if ~isempty(voting_points_start) && ~isempty(voting_points_end)
        flag_do_voting = 1;
    end
end

% Is voting taking place? If so, check those vectors too
if flag_do_voting
    if flag_check_inputs
        % Check the voting_points_start input
        fcn_DebugTools_checkInputsToFunctions(...
            voting_points_start, '2column_of_numbers',Ncircles);
        
        % Check the voting_points_end input
        fcn_DebugTools_checkInputsToFunctions(...
            voting_points_end, '2column_of_numbers',Ncircles);
    end
end



% Does user want to show the plots?
if 8 == nargin
    fig_num = varargin{4};
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
points_innertangent_start = [];
points_innertangent_end = [];
points_outertangent_start = [];
points_outertangent_end = [];


%% Calculate results for inside intersection points
if flag_inside_or_outside<=0
    
    % Calculate inner tangent points
    points_inner = ...
        (centers_end.*radii_start + centers_start.*radii_end)./(radii_start + radii_end);
    
    % Find intersection points using a point and circle formula
    points_innertangent_start = fcn_geometry_findTangentPointsFromPointToCircle(centers_start,radii_start,points_inner);
    points_innertangent_end   = fcn_geometry_findTangentPointsFromPointToCircle(centers_end,radii_end,points_inner);
    
    if flag_do_voting  % This flag is set if user gave voting points
        [points_innertangent_start,start_comparison] = voteGoesToCosest(...
            points_innertangent_start(1:Ncircles,:),...
            points_innertangent_start(Ncircles+1:end,:),...
            voting_points_start);
        [points_innertangent_end,end_comparison] = voteGoesToCosest(...
            points_innertangent_end(1:Ncircles,:),...
            points_innertangent_end(Ncircles+1:end,:),...
            voting_points_end);
        if 1==0
            if ~isequal(start_comparison, end_comparison)
                warning('%s\n','Ambiguous voting detected!',...
                    'The votes for start and end tangents for the inside',...
                    'circles resulted in opposing sides to be selected',...
                    'between the start and ending circles. The result may',...
                    'be outside tangents that do not connect externally.');
                
            end % Ends comparison check
        end
    end % Ends check for voting

    
end % Ends check if inside tangent points should be calculated

%% Calculate outer tangent intersection points. 
if flag_inside_or_outside>=0
    

    % Check to see if these are equi-radii situations. If the test circles
    % have the same radius, then the methods below won't work. So here, we
    % perform special calculations for circles of the same radius, only for
    % data where the radii are equal (noted by the variable equal_indices).
    
    % First, check the radii difference from start and end circle pairs
    radii_difference = radii_start - radii_end;
   
    % flag if radii difference is less than tolerance.
    tol = 0.00001;
    equal_indices = find((radii_difference.^2)<tol^2);
    
    % If any flags exist, do calculations (special) for same circle
    % situations
    if any(equal_indices)
        
        % Calculate the angle of the line connecting the center of the
        % circles
        diff_centers = centers_end-centers_start;
        angles_centers = atan2(diff_centers(:,2),diff_centers(:,1));
        
        % Create 2 new angles, each +/- 90 degrees from that circle-circle
        % angle
        angles_tangents1 = angles_centers+pi/2;
        angles_tangents2 = angles_centers-pi/2;
        
        % Use a polar projection from the center of each circle to
        % calculate where the tangent points would be, using the 2 angles
        % above
        points_outertangent_start_same_radii = [centers_start; centers_start]+...
            [radii_start.*cos(angles_tangents1) radii_start.*sin(angles_tangents1);
             radii_start.*cos(angles_tangents2) radii_start.*sin(angles_tangents2)];
        points_outertangent_end_same_radii = [centers_end; centers_end]+...
            [radii_start.*cos(angles_tangents1) radii_start.*sin(angles_tangents1);
             radii_start.*cos(angles_tangents2) radii_start.*sin(angles_tangents2)];
         
         % Set radii_difference to dummy value to avoid division by zero in
         % the calculations that follow
         radii_difference(equal_indices) = 1;
    end
    
    % Calculate outer tangent points for non-equal-radii case (more typical)
    points_outer = (centers_end.*radii_start - centers_start.*radii_end)./radii_difference;

    % Find intersection points using a point and circle formula (note: this
    % returns 2 times as many points as are in points_outer)
    points_outertangent_start = ...
        fcn_geometry_findTangentPointsFromPointToCircle(...
        centers_start,radii_start,points_outer);
    points_outertangent_end   = ...
        fcn_geometry_findTangentPointsFromPointToCircle(...
        centers_end,radii_end,points_outer);   
    
    % For situations where radii are the same, we used a dummy distance to
    % avoid division by zero, and at the same time calculated the correct
    % answer. Here, we have to substitute the correct answer for these
    % cases back into the solution vector.
    if any(equal_indices) 
        points_outertangent_start(equal_indices,:) = points_outertangent_start_same_radii(equal_indices,:);
        points_outertangent_end(equal_indices,:) = points_outertangent_end_same_radii(equal_indices,:);        
        points_outertangent_start(equal_indices+Ncircles,:) = points_outertangent_start_same_radii(equal_indices+Ncircles,:);
        points_outertangent_end(equal_indices+Ncircles,:) = points_outertangent_end_same_radii(equal_indices+Ncircles,:);        
    end    
    
    % If a vote vector was given, make a vote to keep just one set of
    % tangent points
    if flag_do_voting  % This flag is set if user gave voting points
        [points_outertangent_start,start_comparison] = voteGoesToCosest(...
            points_outertangent_start(1:Ncircles,:),...
            points_outertangent_start(Ncircles+1:end,:),...
            voting_points_start);
        [points_outertangent_end,end_comparison] = voteGoesToCosest(...
            points_outertangent_end(1:Ncircles,:),...
            points_outertangent_end(Ncircles+1:end,:),...
            voting_points_end);
        if ~isequal(start_comparison, end_comparison)
            if 1==0
                warning('%s\n','Ambiguous voting detected!',...
                    'The votes for start and end tangents for the outside',...
                    'circles resulted in opposing sides to be selected',...
                    'between the start and ending circles. The result may',...
                    'be outside tangents that do not connect externally.');
            end
        end % Ends comparison check
    end % Ends check for voting
    
end % Ends check to see if outside tangents should be calculated

% Group inner and outer solutions together, so function can pass these out
points_tangent_start = [points_innertangent_start; points_outertangent_start];
points_tangent_end   = [points_innertangent_end; points_outertangent_end];



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
    fcn_geometry_plotCircle(centers_start,radii_start);
    text(centers_start(:,1),centers_start(:,2)+radii_start*0.5,'S');
    fcn_geometry_plotCircle(centers_end,radii_end);
    text(centers_end(:,1),centers_end(:,2)-radii_start*0.5,'E');
    
    plot(points_tangent_start(:,1),points_tangent_start(:,2),'g+');
    plot(points_tangent_end(:,1),points_tangent_end(:,2),'g+');
    
    if flag_do_voting
        plot(voting_points_start(:,1), voting_points_start(:,2),'co');
        plot(voting_points_end(:,1), voting_points_end(:,2),'co');
    end
    
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

%% For comparisons
function [result,comparison] = voteGoesToCosest(x,y,vote)
% This function takes three vectors of points, x, y, and vote, each is a
% [Nx2] vector of [x y] points that are rowwise separated. The points x or
% y which are closest to vote points, along each row, are returned as the
% result.
do_debug = 0;

num_points = length(x(:,1));
if num_points~=length(y(:,1))
    error('X and Y must be same length');
end
if num_points~=length(vote(:,1))
    error('X and Vote must be same length');
end

if do_debug 
    plot(x(:,1),x(:,2),'ro','Markersize',20);
    plot(y(:,1),y(:,2),'bo','Markersize',20);
    plot(vote(:,1),vote(:,2),'g*','Markersize',20);
end

distances_squared_x_to_voting_point = sum((x - vote).^2,2);
distances_squared_y_to_voting_point = sum((y - vote).^2,2);
comparison = distances_squared_x_to_voting_point < ...
    distances_squared_y_to_voting_point;
not_comparison = ~comparison;
result = x.*comparison + y.*not_comparison;
end

