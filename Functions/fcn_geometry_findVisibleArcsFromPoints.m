function visible_arc_angles = ...
    fcn_geometry_findVisibleArcsFromPoints(...
    centers,...
    radii,...
    points,...
    varargin)
% fcn_geometry_findVisibleArcsFromPoints
% finds the amount of arc visible on a circle from a point external to the
% circle, returning the arc in radians
%
% FORMAT:
%
% visible_arc_angles = ...
%     fcn_geometry_findVisibleArcsFromPoints(...
%     centers,...
%     radii,...
%     points,...
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
%      (OPTIONAL INPUTS)
%
%      fig_num: a figure number to plot results.
%
% OUTPUTS:
%
%      visible_arc_angles: a [N x 1] vector of the arc angle, in radians,
%      for each circle, that is visible from the outside points. NaN is
%      returned for any points that are within a circle.
%
% DEPENDENCIES:
%
%      fcn_geometry_checkInputsToFunctions
%      fcn_geometry_plotCircle
%
% EXAMPLES:
%
% See the script: script_test_fcn_geometry_findVisibleArcsFromPoints
% for a full test suite.
%
% This function was written on 2021_05_22 by S. Brennan
% Questions or comments? sbrennan@psu.edu 

% Revision History:
% 2021-05-23
% -- Created function from fcn_geometry_findTangentPointsFromPointToCircle



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
    narginchk(3,4);
    
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
    
    % Check whether divide by zero will occur, which happens if any points
    % are equal to the center of any circles
    if ~isempty(intersect(points,centers,'rows'))
        error('The query point to calculate a tangent to a circle cannot lie at the center of the same circle');
    end
end


% Does user want to show the plots?
if 4 == nargin
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

%% Preliminary calculations
% Calculate radii squared
radii_squared = radii.^2;

% Calculate distance vector, d, between circle centers and points
d = points - centers;

% Calculate the magnitude of this vector, D, squared
D_squared = sum(d.*d,2);

% Calculate difference in radii
squared_diff = D_squared - radii_squared;
diff = squared_diff.^0.5;

%% Calculate the angles
visible_arc_angles = 2*acos(radii./diff);
    
% For any indices that are bad, set to NaN
visible_arc_angles(squared_diff<0,:) = NaN; 


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
    grid on; grid minor;
        
    % Plot the circles
    fcn_geometry_plotCircle(centers,radii);
    plot(centers(:,1),centers(:,2),'kx');
    axis equal
    
    % Plot the external points
    plot(points(:,1),points(:,2),'go');
    
    % Plot the tangent locations
    points_tangent = ...
        fcn_geometry_findTangentPointsFromPointToCircle(...
        centers,...
        radii,...
        points);
    plot(points_tangent(:,1),points_tangent(:,2),'gx');
    
    % Plot the tangent lines in red, radial lines in magenta
    numpoints = length(points(:,1));
    for i=1:numpoints
        % Tangent lines
        plot([points(i,1) points_tangent(i,1)],[points(i,2) points_tangent(i,2)],'r-');
        plot([points(i,1) points_tangent(i+numpoints,1)],[points(i,2) points_tangent(i+numpoints,2)],'r-');
        
        % Radial spokes
        plot([centers(i,1) points_tangent(i,1)],[centers(i,2) points_tangent(i,2)],'m-');
        plot([centers(i,1) points_tangent(i+numpoints,1)],[centers(i,2) points_tangent(i+numpoints,2)],'m-');

        % Center-to-point
        plot([centers(i,1) points(i,1)],[centers(i,2) points(i,2)],'c--');
        
        % Label the angles
        text(centers(i,1),centers(i,2),sprintf('   visible angle is: %.2f deg',visible_arc_angles(i,1)*180/pi));
    end
end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file); 
end

end % Ends function


