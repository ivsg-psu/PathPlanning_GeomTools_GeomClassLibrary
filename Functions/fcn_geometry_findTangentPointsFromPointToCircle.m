function points_tangent = ...
    fcn_geometry_findTangentPointsFromPointToCircle(...
    centers,...
    radii,...
    points,...
    varargin)
% fcn_geometry_findTangentPointsFromPointToCircle
% finds the two tangent points on a circle given a point.
%
% This calculates the tangent points on circles defined by centers and
% radii, where tangents pass through location given by points. This
% function allows vectorization where centers and points can be
% a vector [x y] where x and y are columns. The radii must be a column
% vector of same length as centers. Points are also in [x y] format, of
% same length as centers. The output

% 
% FORMAT:
%
% points_tangent = ...
%     fcn_geometry_findTangentPointsFromPointToCircle(...
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
%      points_tangent: an [2N x 2] vector of X,Y data for each circle,
%      arranged such that the top points for the N circles appear first,
%      then the bottom points for the N circles. NaN is returned for any
%      points that are within a circle.
%
% DEPENDENCIES:
%
%      fcn_geometry_checkInputsToFunctions
%      fcn_geometry_plotCircle
%
% EXAMPLES:
%
% See the script: script_test_fcn_geometry_findTangentPointsFromPointToCircle
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
% 2021-04-25
% -- Fixed a bug with calculation of intersection
% 2021-05-22
% -- Using external fcn_geometry_plotCircle for plotting



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

% Preliminary calculations
% Calculate radii squared
radii_squared = radii.^2;

% Calculate distance vector, d, between circle centers and points
d = points - centers;

% Calculate the magnitude of this vector, D, squared
D_squared = sum(d.*d,2);

% Calculate difference in radii
squared_diff = D_squared - radii_squared;
diff = squared_diff.^0.5;

% Intermediate calculation requires flipped p2c, 2 versions for 2 points
p2c_flipped1 = [d(:,2) -d(:,1)].*diff;  % The plus-minus solution
p2c_flipped2 = [-d(:,2) d(:,1)].*diff;  % The minus-plus solution

%% Calculate the tangent points
% We do this calculation quickly by stacking the vectors. It's a
% much-abbreviated and faster version than the code listed at the website
% above, but same process.

points_tangent_top = ...
    (radii_squared.*d + radii.*p2c_flipped1)./D_squared + centers;
points_tangent_bottom = ...
    (radii_squared.*d + radii.*p2c_flipped2)./D_squared + centers;

% Check to see if any indices are bad
bad_indices = find(squared_diff<0);
points_tangent_top(bad_indices,:) = NaN; 
points_tangent_bottom(bad_indices,:) = NaN; 

points_tangent = [points_tangent_top; points_tangent_bottom];


% Can do above process all at once via following commands. Same speed and
% less clear.
% points_tangent = [...
%     (radii_squared.*d + radii.*p2c_flipped1);...
%     (radii_squared.*d + radii.*p2c_flipped2)]./[D_squared; D_squared] + ...
%     [centers; centers];
% % Check to see if any indices are bad
% bad_indices = find(~(points_tangent(:,1)==real(points_tangent(:,1))));
% points_tangent(bad_indices,:) = NaN; %#ok<FNDSB>

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
    fcn_geometry_plotCircle(centers,radii);
    plot(centers(:,1),centers(:,2),'kx');
    plot(points(:,1),points(:,2),'go');
    plot(points_tangent(:,1),points_tangent(:,2),'gx');
    numpoints = length(points(:,1));
    for i=1:numpoints
         plot([points(i,1) points_tangent(i,1)],[points(i,2) points_tangent(i,2)],'r-');
         plot([points(i,1) points_tangent(i+numpoints,1)],[points(i,2) points_tangent(i+numpoints,2)],'r-');
    end
end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file); 
end

end % Ends function


