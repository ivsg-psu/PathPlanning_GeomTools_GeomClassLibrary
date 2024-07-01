%% script_test_fcn_geometry_alignArcsInSequence
% Tests the function fcn_geometry_alignArcsInSequence

% Revision history:
% 2024_04_19 - S Brennan
% -- wrote the code
% 2024_05_02 - S Brennan / A. Batchu
% -- add more test cases
% 2024_06_26 - S Brennan 
% -- working with test track data

%% REAL WORLD TEST CASE
% Test with real-world data (test track)



% Check to see if the data was loaded earlier
if ~exist('XY_data','var') 

    % Check to see if XY data for the centerline of the original track lane was loaded earlier
    mat_filename = fullfile(cd,'Data','Centerline_OriginalTrackLane_InnerMarkerClusterCenterOfDoubleYellow.mat');
    if exist(mat_filename,'file')
        load(mat_filename,'XY_data');
    end
end

% Check to see if the fits were calculated earlier
if ~exist('fitSequence_points_forward','var') || ~exist('fitSequence_points_backward','var')

    fitting_tolerance = [30 2]; % Units are meters, in St form

    % Since the XY data is very dense, keep only 1 of every "keep_every" points
    keep_every = 100;
    indicies = (1:length(XY_data(:,1)))';
    small_XY_data_indicies = find(0==mod(indicies,keep_every));
    small_XY_data = XY_data(small_XY_data_indicies,:);

    points = small_XY_data;
    transverse_tolerance = fitting_tolerance(1,2);
    station_tolerance    = fitting_tolerance(1,1);
    points_required_for_agreement = 4;
    flag_force_circle_fit = 0;
    expected_radii_range = [];
    flag_find_only_best_agreement = 0;
    flag_use_permutations = [];

    fig_num = 11;
    figure(fig_num);
    clf;

    domains  = ...
    fcn_geometry_fitHoughCircle(points, transverse_tolerance, ...
            (station_tolerance), (points_required_for_agreement), (flag_force_circle_fit), (expected_radii_range), (flag_find_only_best_agreement), (flag_use_permutations), (fig_num));

    % [fitSequence_points_forward, fitSequence_shapes_forward, fitSequence_endIndicies_forward, fitSequence_parameters_forward, fitSequence_bestFitType_forward] = ...
    %     fcn_geometry_fitSequentialArcs(small_XY_data, fitting_tolerance, flag_fit_backwards, fig_num);

    % Perform the fit forwards and backwards

    flag_fit_backwards = 0;

    fig_num = 12;
    figure(fig_num);
    clf;

    [fitSequence_points_forward, fitSequence_shapes_forward, fitSequence_endIndicies_forward, fitSequence_parameters_forward, fitSequence_bestFitType_forward] = ...
       fcn_geometry_fitSequentialArcsC2(small_XY_data, fitting_tolerance, flag_fit_backwards, fig_num);

    fig_num = 13;
    figure(fig_num);
    clf;

    flag_fit_backwards = 0;
    [fitSequence_points_backward, fitSequence_shapes_backward, fitSequence_endIndicies_backward, fitSequence_parameters_backward, fitSequence_bestFitType_backward] = ...
        fcn_geometry_fitSequentialArcsC2(small_XY_data, fitting_tolerance, 1, fig_num);


end

%% Forward fitting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  ______                               _
% |  ____|                             | |
% | |__ ___  _ ____      ____ _ _ __ __| |
% |  __/ _ \| '__\ \ /\ / / _` | '__/ _` |
% | | | (_) | |   \ V  V / (_| | | | (_| |
% |_|  \___/|_|    \_/\_/ \__,_|_|  \__,_|
%
%
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Forward
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Connect the fits so that the lines perfectly align with the arcs
fig_num = 2;
figure(fig_num);clf;
fitting_tolerance = [10 0.2];
continuity_level = 2;

fits_to_check_types = fitSequence_bestFitType_forward;
fits_to_check_parameters = fitSequence_parameters_forward;

% fits_to_check_types = fitSequence_bestFitType_backward;
% fits_to_check_parameters = fitSequence_parameters_backward;

% clear fits_to_check_types fits_to_check_parameters
% N_fits = 0;
% for ith_fit = [4 5 6]
%     N_fits = N_fits+1;
%     fits_to_check_types{N_fits}      = fitSequence_bestFitType_forward{ith_fit}; %#ok<SAGROW>
%     fits_to_check_parameters{N_fits} = fitSequence_parameters_forward{ith_fit}; %#ok<SAGROW>
% end

[revised_fitSequence_types, revised_fitSequence_parameters, max_feasibility_distance] =  ...
    fcn_geometry_alignArcsInSequence(fits_to_check_types, fits_to_check_parameters, fitting_tolerance, (continuity_level), (fig_num));

test_points_XY = [];
for ith_fit = 1:length(fitSequence_points_forward)
    points_to_plot = fitSequence_points_forward{ith_fit};
    plot(points_to_plot(:,1),points_to_plot(:,2),'k.','MarkerSize',5);
    test_points_XY = [test_points_XY; points_to_plot]; %#ok<AGROW>
end

%% Check the fits
fig_num = 3;
figure(fig_num);
clf;

threshold           = []; %max(max_feasibility_distance,fitting_tolerance(1,2));
curve_test_segment_length = 0.5; % Check every 0.5 meters;

[flag_is_similar, minimum_distance_to_each_point, mean_error, max_error, std_dev_error] = ...
    fcn_geometry_comparePointsToCurve(...
    revised_fitSequence_types, revised_fitSequence_parameters, test_points_XY, ...
    (threshold), (curve_test_segment_length), (fig_num));

%% Backward fitting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  ____             _                           _
% |  _ \           | |                         | |
% | |_) | __ _  ___| | ____      ____ _ _ __ __| |
% |  _ < / _` |/ __| |/ /\ \ /\ / / _` | '__/ _` |
% | |_) | (_| | (__|   <  \ V  V / (_| | | | (_| |
% |____/ \__,_|\___|_|\_\  \_/\_/ \__,_|_|  \__,_|
%
%
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Backward
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Connect the fits so that the lines perfectly align with the arcs
fig_num = 4;
figure(fig_num);clf;
fitting_tolerance = [10 0.2];
continuity_level = 2;

% fits_to_check_types = fitSequence_bestFitType_forward;
% fits_to_check_parameters = fitSequence_parameters_forward;

fits_to_check_types = fitSequence_bestFitType_backward;
fits_to_check_parameters = fitSequence_parameters_backward;

% clear fits_to_check_types fits_to_check_parameters
% N_fits = 0;
% for ith_fit = [4 5 6]
%     N_fits = N_fits+1;
%     fits_to_check_types{N_fits}      = fitSequence_bestFitType_forward{ith_fit}; %#ok<SAGROW>
%     fits_to_check_parameters{N_fits} = fitSequence_parameters_forward{ith_fit}; %#ok<SAGROW>
% end

[revised_fitSequence_types, revised_fitSequence_parameters, max_feasibility_distance] =  ...
    fcn_geometry_alignArcsInSequence(fits_to_check_types, fits_to_check_parameters, fitting_tolerance, (continuity_level), (fig_num));

test_points_XY = [];
for ith_fit = 1:length(fitSequence_points_forward)
    points_to_plot = fitSequence_points_forward{ith_fit};
    plot(points_to_plot(:,1),points_to_plot(:,2),'k.','MarkerSize',5);
    test_points_XY = [test_points_XY; points_to_plot]; %#ok<AGROW>
end

%% Check the fits
fig_num = 5;
figure(fig_num);
clf;

threshold           = []; %max(max_feasibility_distance,fitting_tolerance(1,2));
curve_test_segment_length = 0.5; % Check every 0.5 meters;

[flag_is_similar, minimum_distance_to_each_point, mean_error, max_error, std_dev_error] = ...
    fcn_geometry_comparePointsToCurve(...
    revised_fitSequence_types, revised_fitSequence_parameters, test_points_XY, ...
    (threshold), (curve_test_segment_length), (fig_num));
