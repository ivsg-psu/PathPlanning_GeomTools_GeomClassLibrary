%% script_test_fcn_geometry_alignGeometriesInSequence
% Tests the function fcn_geometry_alignGeometriesInSequence

% Revision history:
% 2024_04_19 - S Brennan
% -- wrote the code
% 2024_05_02 - S Brennan / A. Batchu
% -- add more test cases
% 2024_06_26 - S Brennan 
% -- working with test track data
% 2024_07_21 - S Brennan 
% -- added a few more test cases

% TO-DO
% 2024_07_21 - S. Brennan
% -- Need to add systematic assertions and testing

%% segment joining an offset arc
segment_start = [-20 -0.05];
segment_angle = 0;
segment_length = 60;
segment_parameters = [segment_start segment_angle segment_length];

%               [circleCenter_x.
%                circleCenter_y,
%                radius,
%                start_angle_in_radians, 
%                end_angle_in_radians,
%                flag_this_is_a_circle
%                flag_arc_is_counterclockwise
%               ] 
arc_center = [0 40];
arc_radius = 40;
arc_start_angle = -pi/2;
arc_end_angle = 0;
arc_flag_this_is_a_circle = 0;
arc_flag_arc_is_counterclockwise = 1;
arc_parameters = [arc_center arc_radius arc_start_angle arc_end_angle arc_flag_this_is_a_circle arc_flag_arc_is_counterclockwise];

fits_to_check_types = {'segment','arc'};
fits_to_check_parameters = {segment_parameters, arc_parameters};
fitting_tolerance = 0;
continuity_level = 2;
fig_num = 1;
flag_is_a_loop = 0;

[revised_fitSequence_types, revised_fitSequence_parameters, max_feasibility_distance] =  ...
    fcn_geometry_alignGeometriesInSequence(fits_to_check_types, fits_to_check_parameters, fitting_tolerance, (continuity_level), (flag_is_a_loop), (fig_num));

% Check lengths
assert(length(revised_fitSequence_types)==3);
assert(length(revised_fitSequence_parameters)==3);
assert(length(max_feasibility_distance)==1);

%% REAL WORLD TEST CASE

% Test with real-world data (test track)
fig_num = 1;


% Check to see if the fits were calculated earlier
if ~exist('fitSequence_bestFitType_forward','var') || ~exist('fitSequence_parameters_forward','var')

    % Check to see if XY data for the centerline of the original track lane was loaded earlier
    mat_filename = fullfile(cd,'Data','Centerline_OriginalTrackLane_InnerMarkerClusterCenterOfDoubleYellow.mat');
    if exist(mat_filename,'file')
        load(mat_filename,'XY_data');
    end

    % Since the XY data is very dense, keep only 1 of every "keep_every" points
    keep_every = 20;
    indicies = (1:length(XY_data(:,1)))';
    small_XY_data_indicies = find(0==mod(indicies,keep_every));
    small_XY_data = XY_data(small_XY_data_indicies,:);

    % Perform the fit forwards
    fitting_tolerance = [0.1 10]; % Units are meters
    flag_fit_backwards = 0;
    figure(fig_num);
    clf;

    [fitSequence_points_forward, fitSequence_shapes_forward, fitSequence_endIndicies_forward, fitSequence_parameters_forward, fitSequence_bestFitType_forward] = ...
        fcn_geometry_fitSequentialArcs(small_XY_data, fitting_tolerance, flag_fit_backwards, fig_num);

    % [fitSequence_points_forward, fitSequence_shapes_forward, fitSequence_endIndicies_forward, fitSequence_parameters_forward, fitSequence_bestFitType_forward] = ...
    %    fcn_geometry_fitSequentialArcsC2(small_XY_data, fitting_tolerance, flag_fit_backwards, fig_num);

end


% Connect the fits so that the lines perfectly align with the arcs
fig_num = 2;
figure(fig_num);clf;
fitting_tolerance = [10 0.2];
continuity_level = 2;
flag_is_a_loop = 1;

% fits_to_check_types = fitSequence_bestFitType_forward;
% fits_to_check_parameters = fitSequence_parameters_forward;

clear fits_to_check_types fits_to_check_parameters
N_fits = 0;
for ith_fit = [4 5 6]
    N_fits = N_fits+1;
    fits_to_check_types{N_fits}      = fitSequence_bestFitType_forward{ith_fit}; %#ok<SAGROW>
    fits_to_check_parameters{N_fits} = fitSequence_parameters_forward{ith_fit}; %#ok<SAGROW>
end

[revised_fitSequence_types, revised_fitSequence_parameters, max_feasibility_distance] =  ...
    fcn_geometry_alignGeometriesInSequence(fits_to_check_types, fits_to_check_parameters, fitting_tolerance, (continuity_level), (flag_is_a_loop), (fig_num));

test_points_XY = [];
for ith_fit = 1:length(fitSequence_points_forward)
    points_to_plot = fitSequence_points_forward{ith_fit};
    plot(points_to_plot(:,1),points_to_plot(:,2),'k.','MarkerSize',5);
    test_points_XY = [test_points_XY; points_to_plot]; %#ok<AGROW>
end

% Check the fits
fig_num = 3;
figure(fig_num);
clf;

threshold           = max_feasibility_distance;
curve_test_segment_length = 0.5; % Check every 0.5 meters;

[flag_is_similar, minimum_distance_to_each_point, mean_error, max_error, std_dev_error] = ...
    fcn_geometry_comparePointsToCurve(...
    revised_fitSequence_types, revised_fitSequence_parameters, test_points_XY, ...
    (threshold), (curve_test_segment_length), (fig_num));


