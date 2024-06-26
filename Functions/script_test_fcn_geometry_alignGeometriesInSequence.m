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
    figure(fig_num);
    clf;
    [fitSequence_points_forward, fitSequence_shapes_forward, fitSequence_endIndicies_forward, fitSequence_parameters_forward, fitSequence_bestFitType_forward] = ...
        fcn_geometry_fitSequentialArcsC2(small_XY_data, fitting_tolerance, flag_fit_backwards, fig_num);

end


%% Connect the fits so that the lines perfectly align with the arcs
fig_num = 2;
figure(fig_num);clf;
fitting_tolerance = [100 100];
continuity_level = 2;

clear fits_to_check_types fits_to_check_parameters
% N_fits = 0;
% for ith_fit = [1 2 3 4 5 6]
%     N_fits = N_fits+1;
%     fits_to_check_types{N_fits}      = fitSequence_bestFitType_forward{ith_fit}; %#ok<SAGROW>
%     fits_to_check_parameters{N_fits} = fitSequence_parameters_forward{ith_fit}; %#ok<SAGROW>
% end
fits_to_check_types = fitSequence_bestFitType_forward;
fits_to_check_parameters = fitSequence_parameters_forward;

[revised_fitSequence_types, revised_fitSequence_parameters, max_feasibility_distance] =  ...
    fcn_geometry_alignGeometriesInSequence(fits_to_check_types, fits_to_check_parameters, fitting_tolerance, (continuity_level), (fig_num));

for ith_fit = [1 2 3 5 6]
    points_to_plot = fitSequence_points_forward{ith_fit};
    plot(points_to_plot(:,1),points_to_plot(:,2),'k.','MarkerSize',5);
end
