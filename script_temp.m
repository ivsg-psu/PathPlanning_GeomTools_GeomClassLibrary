%% Now try fitting real-world data
fig_num = 237492;
% Check to see if data was loaded earlier
mat_filename = fullfile(cd,'Data','Centerline_OriginalTrackLane_InnerMarkerClusterCenterOfDoubleYellow.mat');
if exist(mat_filename,'file')
    load(mat_filename,'XY_data');
end


% Perform the fit forwards
fitting_tolerance = 1; % Units are meters
flag_fit_backwards = 0;

% Initialize the subplots
subplot_fig_num = fig_num*200;
trueArcStartIndicies = (1:length(XY_data(:,1)))';
trueNamedCurveTypes = {};
animation_figure_handles = fcn_INTERNAL_setupSubplots(XY_data, trueArcStartIndicies, trueNamedCurveTypes, subplot_fig_num);

[fitSequence_points_forward, fitSequence_shapes_forward, fitSequence_endIndicies_forward, fitSequence_parameters_forward, fitSequence_bestFitType_forward] = ...
    fcn_geometry_fitSequentialArcs(XY_data, fitting_tolerance, flag_fit_backwards, animation_figure_handles, fig_num);

