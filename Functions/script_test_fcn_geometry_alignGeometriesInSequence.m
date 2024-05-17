%% script_test_fcn_geometry_alignGeometriesInSequence
% Tests the function fcn_geometry_alignGeometriesInSequence

% Revision history:
% 2024_04_19 - S Brennan
% -- wrote the code
% 2024_05_02 - S Brennan / A. Batchu
% -- add more test cases


%% REAL WORLD TEST CASE
% Test with real-world data (test track)
fig_num = 1;
figure(fig_num);
clf;

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
    fitting_tolerance = [1 10]; % Units are meters
    flag_fit_backwards = 0;
    [fitSequence_points_forward, fitSequence_shapes_forward, fitSequence_endIndicies_forward, fitSequence_parameters_forward, fitSequence_bestFitType_forward] = ...
        fcn_geometry_fitSequentialArcs(small_XY_data, fitting_tolerance, flag_fit_backwards, 1);

end


%% Connect the fits so that the lines perfectly align with the arcs
fig_num = 2;
figure(fig_num);clf;
fitting_tolerance = [10 2];

clear fits_to_check_types fits_to_check_parameters
for ith_fit = 1:4
    fits_to_check_types{ith_fit}      = fitSequence_bestFitType_forward{ith_fit}; %#ok<SAGROW>
    fits_to_check_parameters{ith_fit} = fitSequence_parameters_forward{ith_fit}; %#ok<SAGROW>
end

fcn_geometry_plotFitSequences(fits_to_check_types, fits_to_check_parameters,(fig_num));

revised_fitSequence_parameters_forward  = ...
    fcn_geometry_alignGeometriesInSequence(fits_to_check_types, fits_to_check_parameters, fitting_tolerance, fig_num);

% fcn_geometry_plotFitSequences(fitSequence_bestFitType_forward, revised_fitSequence_parameters_forward,(fig_num));
