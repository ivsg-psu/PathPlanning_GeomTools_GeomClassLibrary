%% script_test_fcn_geometry_extractXYfromSTSpiral.m
% Exercises the function: fcn_geometry_extractXYfromSTSpiral

% Revision history:
% 2024_03_30 - S. Brennan
% -- wrote the code

close all;



%% BASIC test - perfect fit, perfectly oriented
fig_num = 1;
figure(fig_num); clf;
hold on;

% UNIT CIRCLE ANALYSIS
h0 = 0;
x0 = 0;
y0 = 0;
K0 = 0;

L0 = 2;
s  = (0:0.01:1)'*L0; 
Kf = 1; % Radius of 1

% Call the function fcn_geometry_extractXYfromSTSpiral to predict the
% spiral and calculate the offsets, plotting the results
fcn_geometry_extractXYfromSTSpiral(s,L0,h0,x0,y0,K0,Kf,(fig_num));


% Find the center of the circle tangent at the end of the spiral
% Find the unit vector (need to do this analytically!)
s_tangent = [0.999999 1]'*L0;
[x_tangent,y_tangent] = fcn_geometry_extractXYfromSTSpiral(s_tangent,L0,h0,x0,y0,K0,Kf);
unit_tangent = fcn_geometry_calcUnitVector([diff(x_tangent) diff(y_tangent)]);
unit_orthogonal = unit_tangent*[0 1; -1 0];
calculated_circle_center = (1/Kf)*unit_orthogonal + [x_tangent(end) y_tangent(end)];
y_offset      = calculated_circle_center(1,2) - (1/Kf);

% Plot the circle's center
plot(calculated_circle_center(:,1),calculated_circle_center(:,2),'r+');

% Plot the circle
angles = (0:1:360)'*pi/180;
XY_circle = (1/Kf)*[cos(angles) sin(angles)] + calculated_circle_center;
plot(XY_circle(:,1),XY_circle(:,2),'r-');

% Call the function with no inputs to force a re-calculation of all spiral
% fits. NOTE: NOT YET DONE CODING
% fcn_geometry_findUnitSpiralGivenUnitOffset([],234);

% Now find the offset
[spiralLength, spiralEndAngleInRadians] = fcn_geometry_findUnitSpiralGivenUnitOffset(y_offset, fig_num);

% Check sizes
assert(length(spiralLength(:,1))==1);
assert(length(spiralLength(1,:))==1);
assert(length(spiralEndAngleInRadians(:,1))==1);
assert(length(spiralEndAngleInRadians(1,:))==1);

% Check values
assert(isequal(round(spiralLength,6),round(L0,6)));




%% Test of fast mode
% 
% % Fill circle data
% % circle
% circle_center = [4 3];
% circle_radius = 2;
% M = 3; % 5 points per meter
% sigma = 0.02;
% fig_num = -1;
% 
% circle_test_points = fcn_geometry_fillCircleTestPoints(circle_center, circle_radius, M, sigma); % (fig_num));
% circle_true_parameters = [circle_center, circle_radius, 0, 2*pi, 1];
% 
% % Add outliers?
% % Corrupt the results
% probability_of_corruption = 0.3;
% magnitude_of_corruption = 1;
% 
% corrupted_circle_test_points = fcn_geometry_corruptPointsWithOutliers(circle_test_points,...
%     (probability_of_corruption), (magnitude_of_corruption), (fig_num));
% 
% inputPoints = corrupted_circle_test_points;
% transverse_tolerance = 0.1;
% station_tolerance = [];
% points_required_for_agreement = [];
% flag_force_circle_fit = [];
% expected_radii_range = [];
% flag_find_only_best_agreement = [];
% flag_use_permutations = [];
% 
% 
% domains_corrupted_circle_test_points  = ...
% fcn_geometry_fitHoughCircle(inputPoints, transverse_tolerance, ...
%         (station_tolerance), (points_required_for_agreement), (flag_force_circle_fit), (expected_radii_range), (flag_find_only_best_agreement), (flag_use_permutations), (fig_num));
% 
% % Perform the calculation in slow mode
% fig_num = [];
% REPS = 100; minTimeSlow = Inf;
% tic;
% for i=1:REPS
%     tstart = tic;
% 
%    regression_domain  =  ...
%     fcn_geometry_fitArcRegressionFromHoughFit(domains_corrupted_circle_test_points{1}, fig_num); 
% 
% 
%     telapsed = toc(tstart);
%     minTimeSlow = min(telapsed,minTimeSlow);
% end
% averageTimeSlow = toc/REPS;
% 
% % Perform the operation in fast mode
% fig_num = -1;
% minTimeFast = Inf;
% tic;
% for i=1:REPS
%     tstart = tic;
% 
%     regression_domain  =  ...
%     fcn_geometry_fitArcRegressionFromHoughFit(domains_corrupted_circle_test_points{1}, fig_num); 
% 
%     telapsed = toc(tstart);
%     minTimeFast = min(telapsed,minTimeFast);
% end
% averageTimeFast = toc/REPS;
% 
% fprintf(1,'\n\nComparison of fast and slow modes of fcn_geometry_fitArcRegressionFromHoughFit:\n');
% fprintf(1,'N repetitions: %.0d\n',REPS);
% fprintf(1,'Slow mode average speed per call (seconds): %.8f\n',averageTimeSlow);
% fprintf(1,'Slow mode fastest speed over all calls (seconds): %.8f\n',minTimeSlow);
% fprintf(1,'Fast mode average speed per call (seconds): %.8f\n',averageTimeFast);
% fprintf(1,'Fast mode fastest speed over all calls (seconds): %.8f\n',minTimeFast);
% fprintf(1,'Average ratio of fast mode to slow mode (unitless): %.3f\n',averageTimeSlow/averageTimeFast);
% fprintf(1,'Fastest ratio of fast mode to slow mode (unitless): %.3f\n',minTimeSlow/minTimeFast);


%% Fail conditions
if 1==0
    %% FAIL 1: points not long enough
    points = [2 3];
    [slope,intercept] = fcn_geometry_fitSlopeInterceptNPoints(points,fig_num);
    fprintf(1,'\n\nSlope is: %.2f, Intercept is: %.2f\n',slope,intercept);
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

% function fcn_INTERNAL_printResults(param_string, params)
% header_string = sprintf('%s',param_string);
% fixed_header_string = fcn_DebugTools_debugPrintStringToNCharacters(header_string,25);
% fprintf(1,'%s ',fixed_header_string)
% for ith_value = 1:length(params)
%     param_string = sprintf('%.4f',params(ith_value));
%     fixed_param_string = fcn_DebugTools_debugPrintStringToNCharacters(param_string,15);
%     fprintf(1,'%s ',fixed_param_string)   
% end
% fprintf(1,'\n');
% end