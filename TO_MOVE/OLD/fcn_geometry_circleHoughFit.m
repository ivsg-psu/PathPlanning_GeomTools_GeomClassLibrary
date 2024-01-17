function [fittedParameters, agreementIndices] = fcn_geometry_circleHoughFit(inputPoints, tolerance)
% fcn_geometry_circleHoughFit
%
% This function takes the input points and tolerance as the input and
% outputs the fitted parameters and agreement indices. 
% 
% [fittedParameters, agreementIndices] = fcn_geometry_circleHoughFit(inputPoints, tolerance)
% 
% INPUTS 
% 
% inputPoints: a  Nx2 vector where N is the number of points, but at
%              least 2.
%
% tolerance: The indices points that are in the "tolerance" limit are
% stored. 
%
% OUTPUTS:
%
% fittedParameters: The fitted parameters of the circle are stored here. 
% The size of the matrix is size(nchoosek(1:N,3),1)x2. 
%
% agreementIndices: The indices of the points that are within the tolerance
% limit are stored in this matrix. The size of this matrix is
% size(nchoosek(1:N,3), 1)xN. 
%
% DEPENDENCIES:
%
%      none
%
% EXAMPLES:
%      
% See the script: script_test_fcn_geometry_circleHoughFit
% for a full test suite.
%
% This function was written on 2023_12_15 by Aneesh Batchu
% Questions or comments? abb6486@psu.edu 

% Revision history:
% 2023_12_15 
% -- wrote the code

flag_do_plots = 0; % Flag to plot
flag_do_debug = 0; % Flag to plot the results for debugging
flag_check_inputs = 1; % Flag to perform input checking

if flag_do_debug
    st = dbstack; %#ok<*UNRCH>
    fprintf(1,'STARTING function: %s, in file: %s\n',st(1).name,st(1).file);
    debug_fig_num = 34838;
else
    debug_fig_num = []; %#ok<NASGU>
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
    narginchk(2,2);

end
%% Main Code starts from here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _       
%  |  \/  |     (_)      
%  | \  / | __ _ _ _ __  
%  | |\/| |/ _` | | '_ \ 
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Total number of points
N_points = size(inputPoints,1);

% All possible 3-point combinations
combos_paired = nchoosek(1:N_points,3);
N_combos = size(combos_paired,1);

% Pre-allocation of fittedParameters and agreementIndices for saving
% computation time
fittedParameters = zeros(N_combos,3);
agreementIndices = zeros(N_combos, N_points);

for ith_combo = 1:N_combos
    
    % fcn_INTERNAl_fitPointsTOcircle determines the fitted parameters
    [circleCenter, circleRadius] = fcn_INTERNAl_fitPointsTOcircle(inputPoints(combos_paired(ith_combo,1),:), inputPoints(combos_paired(ith_combo,2),:), inputPoints(combos_paired(ith_combo,3),:));
    % fitted parameters are stored in "fittedParameters" matrix
    fittedParameters(ith_combo,:) = [circleCenter, circleRadius];

    % Finding Agreement Indices

    % Distance of all the input points from the center
    distance_inputpoints_center = sum((inputPoints - circleCenter).^2,2).^0.5;

    % Absolute Error to find the indices in agreement
    abs_error = abs(distance_inputpoints_center - circleRadius); 

    % Indices in agreement
    indicies_in_agreement = (abs_error < tolerance)';
    
    % indices in agreement found in each iteration are stored in
    % "agreementIndices" matrix
    agreementIndices(ith_combo,:) = indicies_in_agreement;

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

end
if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file);
end

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


function [circleCenter, circleRadius] = fcn_INTERNAl_fitPointsTOcircle(point1, point2, point3)

% INPUT
%
% threePoints - x,y coordinates of three points
%
% OUTPUT
%
% circleCenter - center of the fitted circle
% circleRadius - radius of the fitted circle

% For plotting, plots = 1
plots = 0;

% Input must contain three points
threePoints = [point1; point2; point3];

% Displays error if the input does not contain 3 points
if 3 ~= size(threePoints, 1)
    error('Input does not contain three points')
end

% Calculating the center of circles
%
% 1- Find any two chords of the circle - Calculate the distance between two
% point pairs.
% 2- Find the mid-point of the two-point pairs.
% 3- Find the perpendicular bisector of those two chords.
% 4-The intersection of those two perpendicular bisectors is the center of
% the circle.
%
% These formulas are derived based on the geometric properties mentioned
% above.

Denominator = 2 * ( point1(1,1)*(point2(1,2)-point3(1,2)) + point2(1,1)*(point3(1,2)-point1(1,2)) + point3(1,1)*(point1(1,2)-point2(1,2)) );

xcoor_center = ( (point1(1,1)^2+point1(1,2)^2)*(point2(1,2)-point3(1,2)) + ...
    (point2(1,1)^2+point2(1,2)^2)*(point3(1,2)-point1(1,2)) + ...
    (point3(1,1)^2+point3(1,2)^2)*(point1(1,2)-point2(1,2)) ) / Denominator;

ycoor_center = ( (point1(1,1)^2+point1(1,2)^2)*(point3(1,1)-point2(1,1)) + ...
    (point2(1,1)^2+point2(1,2)^2)*(point1(1,1)-point3(1,1)) + ...
    (point3(1,1)^2+point3(1,2)^2)*(point2(1,1)-point1(1,1)) ) / Denominator;

circleCenter = [xcoor_center, ycoor_center];

% Circle radius is obtained by calculating the distance between the center
% and any point on the circle.

circleRadius = ( (point1(1,1)-xcoor_center)^2 + (point1(1,2)-ycoor_center)^2 )^0.5;

if plots

    % Plot the points and the fitted circle
    % theta_start = atan2(point1(1,2) - circleCenter(1,2), point1(1,1) - circleCenter(1,1));
    % theta_end = atan2(point3(1,2) - circleCenter(1,2), point3(1,1) - circleCenter(1,1));
    % theta = linspace(theta_start, theta_end, 100);

    theta = linspace(0, 2*pi, 100);
    circle_xcoor = circleCenter(1,1) + circleRadius*cos(theta);
    circle_ycoor = circleCenter(1,2) + circleRadius*sin(theta);

    % plot the circle
    figure(10023)
    plot(circle_xcoor, circle_ycoor, 'r', 'LineWidth',2);
    hold on

    % Plot the threePoints on the circle
    scatter(threePoints(:,1), threePoints(:,2), 'b','filled');

    % Plot the center
    scatter(circleCenter(1,1), circleCenter(1,2), 'MarkerEdgeColor', 'g', 'MarkerFaceColor', 'g');

    % Details
    xlabel('X-axis');
    ylabel('Y-axis');
    title('Fitted Circle Passing Through Three Points');
    legend('Fitted Circle', 'Points', 'Center');

    axis equal
    box on
    grid on
    
    % texts = 0, for no text on figure
    texts = 1;

    if texts
        % Three points
        text(point1(1,1), point1(1,2), '  Point 1');
        text(point2(1,1), point2(1,2), '  Point 2');
        text(point3(1,1), point3(1,2), '  Point 3');
        % Center
        text(circleCenter(1,1), circleCenter(1,2), ['  Center (' num2str(circleCenter(1,1)) ', ' num2str(circleCenter(1,2)) ')' ])
    end

end

end