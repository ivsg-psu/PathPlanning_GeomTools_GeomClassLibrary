%% script_test_fcn_geometry_isFeasibleAlignGeomSeries
% Tests the function fcn_geometry_isFeasibleAlignGeomSeries

% Revision history:
% 2024_05_18 - S Brennan
% -- wrote the code
close all;

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
    fitting_tolerance = [0.5 10]; % Units are meters
    flag_fit_backwards = 0;
    
    figure(fig_num);
    plot(XY_data(:,1),XY_data(:,2),'k.','MarkerSize',20)
    
    [fitSequence_points_forward, fitSequence_shapes_forward, fitSequence_endIndicies_forward, fitSequence_parameters_forward, fitSequence_bestFitType_forward] = ...
        fcn_geometry_fitSequentialArcs(small_XY_data, fitting_tolerance, flag_fit_backwards, fig_num);


end


%% Check if the fits are feasible
fig_num = 2;
figure(fig_num);clf;
fitting_tolerance = 0;

clear fits_to_check_types fits_to_check_parameters
fits_to_check_types = fitSequence_bestFitType_forward;
fits_to_check_parameters = fitSequence_parameters_forward;
% for ith_fit = 1:6
%     fits_to_check_types{ith_fit}      = fitSequence_bestFitType_forward{ith_fit}; %#ok<SAGROW>
%     fits_to_check_parameters{ith_fit} = fitSequence_parameters_forward{ith_fit}; %#ok<SAGROW>
% end

plot(XY_data(:,1),XY_data(:,2),'k.','MarkerSize',20)
fcn_geometry_plotFitSequences(fits_to_check_types, fits_to_check_parameters,(fig_num));

distances_by_continuity = zeros(length(fits_to_check_types)-1,3);
for ith_continuity = 1:3
    continuity_level = ith_continuity-1;

    [joint_feasibility, joint_distance_of_failure] = ...
        fcn_geometry_isFeasibleAlignGeomSeries(fits_to_check_types, fits_to_check_parameters, continuity_level,(fitting_tolerance),(fig_num));

    distances_by_continuity(:,ith_continuity) = joint_distance_of_failure;

end
disp(distances_by_continuity);