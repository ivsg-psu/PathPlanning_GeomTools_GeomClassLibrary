%% script_test_fcn_geometry_curvatureAlongCurve
% Exercises the function: fcn_geometry_curvatureAlongCurve

% 2024_06_27 - S. Brennan
% -- wrote the code

close all;

flag_be_verbose = 1;

%% Fill Data
rng(1); % Fix the random number, for debugging

% Set 1==1 to use test track data. Set 1==0 to use artificial data.
if 1==1
    % TEST TRACK DATA

    % Check to see if XY data for the centerline of the original track lane was loaded earlier
    mat_filename = fullfile(cd,'Data','Centerline_OriginalTrackLane_InnerMarkerClusterCenterOfDoubleYellow.mat');
    if exist(mat_filename,'file')
        load(mat_filename,'XY_data');
    end

    % Pre-append and post-append data, to allow wrap-around?
    % ADD THIS?

    % Since the XY data is very dense, keep only 1 of every "keep_every" points
    keep_every = 20; % 100 works OK
    indicies = (1:length(XY_data(:,1)))';
    small_XY_data_indicies = find(0==mod(indicies,keep_every));
    small_XY_data = XY_data(small_XY_data_indicies,:);
    points_to_fit = small_XY_data;
else
    % ARTIFICIAL ARC DATA

    % arc_pattern has [1/R and L] for each segment as a row
    arc_pattern = [...
        1/20, 15;
        0 20;
        -1/5 10;
        0 10;
        1/15 40;
        0 15
        -1/10 20];

    % arc_pattern = [...
    %     1/20, 45;
    %     0 20;
    %     -1/10 10;
    %     0 10];

    M = 5; % How many points per meter
    sigma = 0.02; % The standard deviation in the points relative to the perfect function fit, in meters

    [points_to_fit, ~, ~, trueArcStartIndicies, trueNamedCurveTypes, trueParameters] = fcn_geometry_fillArcSequenceTestPoints(arc_pattern, M, sigma, -1);

end

% Add the tiniest bit of noise - this prevents singularities when doing
% regressions
points_to_fit = points_to_fit + 0.001*rand(size(points_to_fit));

% Plot the raw data
fig_num = 1;
figure(fig_num);
clf;

hold on;
grid on;
axis equal
xlabel('X [m]');
ylabel('Y [m]');

plot(points_to_fit(:,1),points_to_fit(:,2),'k.','MarkerSize',20);

%% Assess the data properties
% to find average spacing, number of points, etc.
Npoints = length(points_to_fit(:,1));
minimum_island_separation = 10; % Units are meters. See comments below for explanation

spatial_differences = diff(points_to_fit(:,1:2));
spatial_distances   = real(sum(spatial_differences.^2,2).^0.5);
average_spacing     = mean(spatial_distances);

% Check if the data is a loop (this changes a few steps that follow)
distance_XY_from_start = real(sum((points_to_fit(:,1:2) - points_to_fit(1,1:2)).^2,2).^0.5);
max_distance_XY_from_start = max(distance_XY_from_start);
distance_start_to_end  = real(sum((points_to_fit(end,1:2) - points_to_fit(1,1:2)).^2,2).^0.5);
flag_is_a_loop = 0;

% Make sure the start/end are not 5 standard deviations away from each
% other, and no more than 2 island separations (e.g. there's a chance they
% are connected to each other).
if distance_start_to_end<(5*average_spacing) || distance_start_to_end<(2*minimum_island_separation)
    flag_is_a_loop = 1;
end

if 1==flag_be_verbose
    % Report results
    fprintf(1,'\n\n');
    fprintf(1,'RAW DATA ASSESSMENT:\n');
    fprintf(1,'Number of points: %.0f\n',Npoints);
    fprintf(1,'Spacing:\n');
    fprintf(1,'Mean spacing: %.3f meters\n',average_spacing);
    fprintf(1,'Stdv spacing: %.3f meters\n',std(spatial_distances));
    fprintf(1,'Max  spacing: %.3f meters\n',max(spatial_distances));
    fprintf(1,'Min  spacing: %.3f meters\n',min(spatial_distances));
    fprintf(1,'Looping:\n');
    fprintf(1,'Max data distance from start: %.3f meters\n',max_distance_XY_from_start);
    fprintf(1,'Distance from start to end: %.3f meters\n',distance_start_to_end);
    fprintf(1,'Is this data a loop? %.0f \n',flag_is_a_loop);
end



%% Show how to use fcn_geometry_curvatureAlongCurve to find any islands in geometric information
% Islands are where there are interconnected arcs that have no line
% segments within. Islands are useful for analysis because calculations
% done within an island are not affected by calculations in other islands,
% and so the data can be sub-grouped by island and processed island by
% island, thereby saving huge amounts of computation.
%
% We find islands by calculating if there are "gaps" in the
% data wherein the curvature of the data, within the gaps, is
% indistinguishable from a line fit. This calculation is within the
% curvatureAtPoint function. The curvature calculation is VERY slow as it
% checks, via regression, all possible circles that can be created with the
% given point at the center of the data set. One can set number of data
% points to consider to right and left of the test point. Default is to use
% them all, but for large data this is VERY slow. Instead of using them
% all, we can specify the minimum distance ever expected between islands,
% e.g. the minimum distance allowed to be a "straightaway" on a road. This
% is a user-defined parameter that of course depends on the actual road.
% From this, we can calculate how many data points should be used for
% curvature calculations. Of note: data_width needs to be at least 2 or
% more for the curvature calculation to work.

% Set up the distances wherein the islands must be separated to be an
% "island"
if Npoints>100
    data_width = ceil(minimum_island_separation/average_spacing);
    if data_width<=2
        data_width = [];
        warning('Not enough data for curvature calculations. Using all data.');
    end

    % Enforce a minimum number of points, otherwise the SNR is very poor- we do not want to use less than
    % 20 points typically
    data_width = max(data_width,20);
end

fig_num = 1234;


[curvatures, arc_centers, index_ranges, point_curvature_minimum, curvature_SNRs] = fcn_geometry_curvatureAlongCurve(points_to_fit, (data_width), (fig_num));

% Assign islands to locations where the SNR is less than 1, e.g. it's more
% likely that the data is a line than an arc. We make it 3 here to give it
% a bit of wiggle-room (some are right on edge).
is_island = curvature_SNRs>3;

%% Show how to use fcn_geometry_curvatureAlongCurve to fix loops
% In the test track, and in any data that forms a loop, it is common that
% the start/end of data will be on a curve. The results will be an island
% at start and at end, but with straightaways somehwere before and after
% the end. 
%
% If this is the case, they can be fixed by shifting the data to start at
% the middle of the first "lake" part after the starting island. This is
% basically "fixing" the loops so that the loop's start/end does not break
% a continuous island of arcs in half by putting the start/end point in the
% middle of an arc (as is the case at the test track).

lake_exists_in_data = any(~is_island);
if is_island(1,1) && is_island(end,1) && lake_exists_in_data
    % Need to fix the data

    island_starts = find(~is_island,1);
    % Make everything a lake up to where the lake starts. This makes
    % finding where the lake ends very trivial
    filled_data = is_island;
    filled_data(1:island_starts) = 0; 
    lake_ends = find(filled_data==1,1);
    mid_lake = round((island_starts+lake_ends)/2);

    % Shift the points so they start mid-lake
    shifted_points = [points_to_fit(mid_lake+1:end,:); points_to_fit(1:mid_lake,:)];
    
    % Recalculate the data using the shifted points
    fig_num = 2345;
    [curvatures, arc_centers, index_ranges, point_curvature_minimum, curvature_SNRs] = fcn_geometry_curvatureAlongCurve(shifted_points, (data_width), (fig_num));
else
    shifted_points = points_to_fit;
end
is_island = curvature_SNRs>3;

% Plot the results
fit_fig_num = 13232;
figure(fit_fig_num);
clf;

hold on;
grid on;
axis equal
xlabel('X [m]');
ylabel('Y [m]');

plot(shifted_points(:,1),shifted_points(:,2),'k.','MarkerSize',20);
plot(shifted_points(is_island,1),shifted_points(is_island,2),'r.','MarkerSize',20);

%% Use curvatureGroupAssignment to separate contiguous points into islands
fig_num = 67543;
island_ranges = fcn_INTERNAL_curvatureGroupAssignment(shifted_points, is_island, fig_num);

%% Use extractModelsFromCurvature to convert each island into arcs
entire_model_fit_ID_at_each_index = nan(Npoints,1);
entire_arc_matrix_C2 = [];
entire_model_SNRs = [];
Nmodels = 0;

SNR_threshold = 30;
for ith_island = 1:length(island_ranges)
    this_island_range = island_ranges{ith_island};
    this_island_points = shifted_points(this_island_range,:);

    if 1==0
        [~, ~, model_SNRs, sorted_model_fit_ID_at_each_index] = fcn_INTERNAL_extractModelsFromCurvature(this_island_points, SNR_threshold);
    else
        % STEP1: Calculate the full curvatures and SNRs of this data
        debug_fig_num = 2346;
        data_width = []; % Default is to use all possible points
        [curvatures, arc_centers, index_ranges, point_curvature_minimum, curvature_SNRs] = fcn_geometry_curvatureAlongCurve(this_island_points, (data_width), (debug_fig_num));

        % STEP2: Use the curvature SNRs to extract models
        fig_num = 75655;
        [best_fit_arcs, best_fit_SNRs, best_fit_ranges, model_fit_ID_at_each_index] = fcn_INTERNAL_curvatureArcsFromSNR(this_island_points, curvature_SNRs, SNR_threshold, arc_centers, curvatures, point_curvature_minimum, index_ranges, fig_num);

        % STEP3: order the models to be in correct order
        fig_num = 383834;
        all_segments = [];
        [~, ~, model_SNRs, sorted_model_fit_ID_at_each_index, arc_matrix, ~] = fcn_INTERNAL_curvatureModelSort(best_fit_arcs, all_segments, best_fit_SNRs, model_fit_ID_at_each_index, fig_num);

        % STEP4: ensure arcs have C2 continuity
        fig_num = 75633;
        arc_matrix_C2 = fcn_INTERNAL_alignArcsBySNRandC2(arc_matrix, model_SNRs, sorted_model_fit_ID_at_each_index, fig_num);

    end

    % Update the model IDs
    % Set any unfilled (zero) values to NaN
    sorted_model_fit_ID_at_each_index(sorted_model_fit_ID_at_each_index==0) = nan;

    % Offset all the model IDs that were just measured so they match the
    % rows of the updated parameter list
    offset_model_fit_ID_at_each_index = Nmodels + sorted_model_fit_ID_at_each_index;

    % Copy these data into the correct range area.
    entire_model_fit_ID_at_each_index(this_island_range) = offset_model_fit_ID_at_each_index;

    % Update the parameter lists
    entire_arc_matrix_C2 = [entire_arc_matrix_C2; arc_matrix_C2]; %#ok<AGROW>
    if ~isempty(entire_arc_matrix_C2)
        Nmodels = length(entire_arc_matrix_C2(:,1));
    end

    % Update the SNRs
    entire_model_SNRs = [entire_model_SNRs, model_SNRs];

end

%% Fill the unfilled areas with line segments
NumFitsGood = length(entire_arc_matrix_C2(:,1));
revised_entire_model_fit_ID_at_each_index = entire_model_fit_ID_at_each_index;

indicies_unfilled = isnan(entire_model_fit_ID_at_each_index);

% Initialize the segments
all_segments = [];

flag_need_to_fill_start = 0;

Nsegments = 0;
while any(indicies_unfilled)

    start_index = find(indicies_unfilled==1,1);
    remainder   = indicies_unfilled;
    remainder(1:start_index) = 1; % Fill in any zeros beforehand
    end_index   = find(remainder==0,1); % Find first zero afterwards



    % Check for situation where a loop is causing the first and last areas
    % to be segments
    flag_skip_fit = 0;
    if 1==start_index && flag_is_a_loop
        % Grab the arc at end, and use it as the start
        start_index = find(~indicies_unfilled==1,1,'last');
        fit_before = entire_model_fit_ID_at_each_index(start_index);
        fit_after  = entire_model_fit_ID_at_each_index(end_index);

        % Shut off the indicies searched thus far, in prep for next round,
        % including the wrap-around portion
        indicies_unfilled(1:end_index) = 0;
        indicies_unfilled(start_index:end) = 0;

    elseif 1==start_index
        % Shut off the indicies searched thus far, in prep for next round
        indicies_unfilled(1:end_index) = 0;
        flag_skip_fit = 1;
    else
        % Grab the arc before and after the open area
        fit_before = entire_model_fit_ID_at_each_index(start_index-1);
        fit_after  = entire_model_fit_ID_at_each_index(end_index);

        % Shut off the indicies searched thus far, in prep for next round
        indicies_unfilled(1:end_index) = 0;

    end


    % Check for errors
    if isempty(end_index)
        error('Empty end indicies should never occur')
    end

    if 0==flag_skip_fit
        Nsegments = Nsegments+1;

        arc_before = entire_arc_matrix_C2(fit_before,:);
        arc_after  = entire_arc_matrix_C2(fit_after,:);

        continuity_level = 1;
        [revised_arc1_parameters, revised_arc2_parameters, revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters]  = ...
            fcn_geometry_alignArcArc(arc_before, arc_after, (transverse_tolerance), (continuity_level),  (57894));

        % Save arc revisions
        entire_arc_matrix_C2(fit_before,:) = revised_arc1_parameters; %#ok<SAGROW>
        entire_arc_matrix_C2(fit_after,:) = revised_arc2_parameters; %#ok<SAGROW>

        segment_start_XY = revised_intermediate_geometry_join_parameters(1,1:2);
        segment_angle    = revised_intermediate_geometry_join_parameters(1,3);
        segment_length   = revised_intermediate_geometry_join_parameters(1,4);
        segment_end_XY   = segment_start_XY + [cos(segment_angle) sin(segment_angle)]*segment_length;

        % Make sure a segment is produced
        assert(strcmp(revised_intermediate_geometry_join_type,'segment'));

        % Shift the segment away from the arcs very slightly so that C2
        % continuity can be recovered
        arc_before_center = arc_before(1,1:2);
        arc_after_center  = arc_after(1,1:2);

        vector_from_arc_before_center_to_segment_start = segment_start_XY - arc_before_center;
        vector_from_arc_after_center_to_segment_end    = segment_end_XY   - arc_after_center;
        unit_vector_from_arc_before_center_to_segment_start = fcn_geometry_calcUnitVector(vector_from_arc_before_center_to_segment_start);
        unit_vector_from_arc_after_center_to_segment_end    = fcn_geometry_calcUnitVector(vector_from_arc_after_center_to_segment_end);

        offset_nudge = 0.05;
        revised_segment_start_XY = segment_start_XY + unit_vector_from_arc_before_center_to_segment_start*offset_nudge;
        revised_segment_end_XY   = segment_end_XY   + unit_vector_from_arc_after_center_to_segment_end*offset_nudge;
        revised_segment_vector   = revised_segment_end_XY - revised_segment_start_XY;
        revised_segment_angle    = atan2(revised_segment_vector(2),revised_segment_vector(1));
        revised_segment_length   = real(sum(revised_segment_vector.^2,2).^0.5);
        revised_segment          = [revised_segment_start_XY revised_segment_angle revised_segment_length];

        % Save results
        all_segments = [all_segments; revised_segment]; %#ok<AGROW>

        % Indicate where the line segments (negative numbering) are:
        if start_index>end_index && flag_is_a_loop
            revised_entire_model_fit_ID_at_each_index(start_index+1:end) = -1*Nsegments;
            revised_entire_model_fit_ID_at_each_index(1:end_index-1) = -1*Nsegments;

        else
            revised_entire_model_fit_ID_at_each_index(start_index:(end_index-1),1) = -1*Nsegments;
        end
    end
end

% Plot the revised results
arc_line_fig_num = 23434;
figure(arc_line_fig_num);
clf;

hold on;
grid on;
axis equal
xlabel('X [m]');
ylabel('Y [m]');

for ith_fit = 1:length(entire_arc_matrix_C2(:,1))
    fcn_geometry_plotGeometry('arc',entire_arc_matrix_C2(ith_fit,:));
end


for ith_fit = 1:length(all_segments(:,1))
    fcn_geometry_plotGeometry('segment',all_segments(ith_fit,:));
end

%% Sort the fits to be in the correct order
fig_num = 1111;
[fitSequence_fitTypes, fitSequence_parameters, model_SNRs, sorted_model_fit_ID_at_each_index, arc_matrix, segment_matrix] = fcn_INTERNAL_curvatureModelSort(entire_arc_matrix_C2, all_segments, [], revised_entire_model_fit_ID_at_each_index, fig_num);


%% Align fits to each other to ensure C2 continuity
fig_num = 76767;
continuity_level = 2;
[revised_fitSequence_types, revised_fitSequence_parameters] = ...
    fcn_geometry_alignGeometriesInSequence(fitSequence_fitTypes, fitSequence_parameters, transverse_tolerance, (continuity_level), (fig_num));

%% Analyze error of model fitting
data_to_analyze = points_to_fit; % small_XY_data
plot(data_to_analyze(:,1),data_to_analyze(:,2),'k.','MarkerSize',5);

% Check the fits
fig_num = 1000*ith_fit;
figure(fig_num);
clf;

threshold           = []; %max(max_feasibility_distance,fitting_tolerance(1,2));
curve_test_segment_length = 0.5; % Check every 0.5 meters;

% [flag_is_similar, minimum_distance_to_each_point, mean_error_forward, max_error_forward, std_dev_error_forward] = ...
[flag_is_similar, minimum_distance_to_each_point, indicies_of_nearest_reference_points, mean_error, max_error, std_dev_error] = ...
    fcn_geometry_comparePointsToCurve(...
    revised_fitSequence_types, revised_fitSequence_parameters, data_to_analyze, ...
    (threshold), (curve_test_segment_length), (fig_num));





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

%% fcn_INTERNAL_curvatureGroupAssignment
function island_ranges = fcn_INTERNAL_curvatureGroupAssignment(shifted_points, is_island, fig_num)
% Given an [Nx2] array of input points and an [Nx1] array flagging each
% point as an island or not, this function groups the contigous groups or
% "islands" of points into separate ranges, returning the range for each
% "island" grouping

% Find the number of islands
jumps_in_islands = [is_island(1,1); diff(is_island); 0];

% Each transition "up" from 0 to 1 will start a new island. We can count
% these to determine number of island starts. Similarly, each transition
% "down" from 1 to 0 will end the island. Again, we can count these

indicies_island_start = find(jumps_in_islands>0.5);
indicies_island_end = find(jumps_in_islands<-0.5);

% Make sure these match in length!
assert(length(indicies_island_start)==length(indicies_island_end));

% How many islands do we have in this data?
N_islands = length(indicies_island_start);

% Initialize the output
island_ranges{N_islands} = [];

% Fill the outputs
for ith_island = 1:N_islands
    island_ranges{ith_island} = (indicies_island_start(ith_island):indicies_island_end(ith_island)); 
end

% Plot the results?
if ~isempty(fig_num)
    figure(fig_num);
    clf;

    hold on;
    grid on;
    axis equal
    xlabel('X [m]');
    ylabel('Y [m]');

    % Plot the points
    plot(shifted_points(:,1),shifted_points(:,2),'k.','MarkerSize',20);

    % Plot the islands
    for ith_island = 1:N_islands
        this_island_range = island_ranges{ith_island};
        plot(shifted_points(this_island_range,1),shifted_points(this_island_range,2),'.','MarkerSize',10);
    end
end
end % Ends fcn_INTERNAL_curvatureGroupAssignment

%% fcn_INTERNAL_extractModelsFromCurvature
function [fitSequence_fitTypes, fitSequence_parameters, model_SNRs, sorted_model_fit_ID_at_each_index] = fcn_INTERNAL_extractModelsFromCurvature(this_island_points, SNR_threshold)
% Steps:
% For each island, extractModelsFromCurvature function does the following
% STEP1:  find full curvatures and SNRs
% STEP2:  use the curvature SNRs to extract models at each island, recording model at each index
% STEP3:  order the models in each island to be in correct order

% STEP4?      correct the models in each island to ensure C2 continuity
%% STEP1: Calculate the full curvatures and SNRs of this data
debug_fig_num = 2346;
data_width = []; % Default is to use all possible points
[curvatures, arc_centers, index_ranges, point_curvature_minimum, curvature_SNRs] = fcn_geometry_curvatureAlongCurve(this_island_points, (data_width), (debug_fig_num));

%% STEP2: Use the curvature SNRs to extract models
fig_num = 75655;
[best_fit_arcs, best_fit_SNRs, best_fit_ranges, model_fit_ID_at_each_index] = fcn_INTERNAL_curvatureArcsFromSNR(this_island_points, curvature_SNRs, SNR_threshold, arc_centers, curvatures, point_curvature_minimum, index_ranges, fig_num);

%% STEP3: order the models to be in correct order
fig_num = 383834;
all_segments = [];
[fitSequence_fitTypes, fitSequence_parameters, model_SNRs, sorted_model_fit_ID_at_each_index, arc_matrix, segment_matrix] = fcn_INTERNAL_curvatureModelSort(best_fit_arcs, all_segments, best_fit_SNRs, model_fit_ID_at_each_index, fig_num);

%% STEP4: ensure arcs have C2 continuity
fig_num = 75633;
[fitSequence_fitTypes_C2, fitSequence_parameters_C2] = fcn_INTERNAL_alignArcsBySNRandC2(fitSequence_fitTypes, fitSequence_parameters, model_SNRs, sorted_model_fit_ID_at_each_index, fig_num);

end % Ends fcn_INTERNAL_extractModelsFromCurvature


%% fcn_INTERNAL_curvatureArcsFromSNR
function [best_fit_arcs, best_fit_SNRs, best_fit_ranges, model_fit_ID_at_each_index] = fcn_INTERNAL_curvatureArcsFromSNR(this_island_points, curvature_SNRs, SNR_threshold, arc_centers, curvatures, point_curvature_minimum, index_ranges, fig_num)

% Find number of points
Npoints = length(this_island_points(:,1));
indiciesPoints = (1:Npoints)';

% Initialize output variables
best_fit_arcs = [];
best_fit_SNRs = [];
best_fit_ranges = [];

% Initialize internal array that keeps track of which SNRs we have already
% analyzed. If this array has an NaN value, then that point has been
% analyzed or should be ignored.
remaining_curvature_SNRs = curvature_SNRs;
remaining_curvatures = curvatures;

% Ignore any points that have an SNR below the threshold
remaining_curvature_SNRs(curvature_SNRs<SNR_threshold) = nan;
remaining_curvatures(curvature_SNRs<SNR_threshold) = nan;

% Keep track of which models go with which points by setting all the model
% numbers to zero
model_fit_ID_at_each_index = 0*curvature_SNRs;

% Set up a debugging figure
flag_do_debug = 1;
debug_fig_num = 38383;
figure(debug_fig_num);clf;

flag_first_time = 1;

NumFits = 0;

while ~all(isnan(remaining_curvature_SNRs))

    % Find the best remaining
    [~,best_SNR_index] = max(remaining_curvature_SNRs);
    NumFits = NumFits+1;
    
    % Find the index range of this fit
    min_index = best_SNR_index-index_ranges(best_SNR_index);
    max_index = best_SNR_index+index_ranges(best_SNR_index);
    this_index_range = (min_index:max_index)';

    % Plot results?
    if flag_do_debug
        %%%%%%%%%%%%%%%%%%%%%%
        figure(debug_fig_num)
        subplot(1,3,1);
        % cla;

        hold on;
        grid on;
        axis equal
        xlabel('X [m]');
        ylabel('Y [m]');

        % Make axis slightly larger?
        if 1==flag_first_time

            % Plot the input points
            plot(this_island_points(:,1),this_island_points(:,2),'b.','MarkerSize',20);


            temp = axis;
            axis_range_x = temp(2)-temp(1);
            axis_range_y = temp(4)-temp(3);
            percent_larger = 0.3;
            axis([temp(1)-percent_larger*axis_range_x, temp(2)+percent_larger*axis_range_x,  temp(3)-percent_larger*axis_range_y, temp(4)+percent_larger*axis_range_y]);
            temp_axis = axis;
            flag_first_time = 0;
        else
            axis(temp_axis);
        end

        % Plot the circle fit at the point
        fcn_geometry_plotCircle(arc_centers(best_SNR_index,:), 1/curvatures(best_SNR_index),'r-',debug_fig_num);


        % Plot the index range
        plot(this_island_points(min_index:max_index,1),this_island_points(min_index:max_index,2),'m.','MarkerSize',10)

        % Plot the max SNR point
        plot(this_island_points(best_SNR_index,1),this_island_points(best_SNR_index,2),'g.','MarkerSize',30)

        axis(temp_axis);

        title('Input points');

        %%%%%%%%%%%%%%%%%%%%%%
        subplot(1,3,2);
        % cla;

        semilogy(indiciesPoints,curvatures,'k-');
        hold on;
        semilogy(indiciesPoints,remaining_curvatures,'b-', 'LineWidth',2);
        semilogy(indiciesPoints,point_curvature_minimum,'-','Color',[0.6 0.6 0.6]);

        % Plot the part covered by this fit
        plot(indiciesPoints(this_index_range),curvatures(best_SNR_index)*ones(length(this_index_range)),'k-','Markersize',20,'LineWidth',3);

        grid on;
        xlabel('Index [count]');
        ylabel('Curvature [1/m]');
        title('Curvatures')

        %%%%%%%%%%%%%%%%%%%%%%
        subplot(1,3,3);
        % cla;
        grid on;
        hold on;

        plot(indiciesPoints, curvature_SNRs,'k-');
        plot(indiciesPoints, remaining_curvature_SNRs,'b-', 'LineWidth',2);

        % Plot the part covered by this fit
        plot(indiciesPoints(this_index_range), remaining_curvature_SNRs(this_index_range,:),'k-','LineWidth',3);

        % Plot the best point
        plot(indiciesPoints(best_SNR_index), remaining_curvature_SNRs(best_SNR_index),'g.','Markersize',30);

        xlabel('Index [count]');
        ylabel('SNR [unitless]');
        title('Curvature SNR')
    end
    
    % Save results
    best_fit_SNRs    = [best_fit_SNRs; remaining_curvature_SNRs(best_SNR_index)]; %#ok<AGROW>
    best_fit_ranges = [best_fit_ranges; index_ranges(best_SNR_index)]; %#ok<AGROW>
    

    % Block out the indicies of this fit
    remaining_curvature_SNRs(this_index_range) = nan;

    % Check which portion of these points has not yet been filled (model is
    % 0) and ONLY fill these with the current model fit number
    is_zeros = model_fit_ID_at_each_index(this_index_range)==0;
    unfilled_index_range = this_index_range(is_zeros);
    model_fit_ID_at_each_index(unfilled_index_range) = NumFits;

    % Save fit parameters. 
    % The standard arc parameter format:
    %               [circleCenter_x.
    %                circleCenter_y,
    %                radius,
    %                start_angle_in_radians,
    %                end_angle_in_radians,
    %                flag_this_is_a_circle
    %                flag_arc_is_counterclockwise
    %               ]
    best_arc_center                    = arc_centers(best_SNR_index,:);
    best_arc_radius                    = 1/curvatures(best_SNR_index);
    vector_from_circle_center_to_start = this_island_points(min_index,:)-best_arc_center;
    vector_from_circle_center_to_end   = this_island_points(max_index,:)-best_arc_center;
    best_arc_start_angle_in_radians    = mod(atan2(vector_from_circle_center_to_start(2),vector_from_circle_center_to_start(1)),2*pi);
    best_arc_end_angle_in_radians      = mod(atan2(vector_from_circle_center_to_end(2),vector_from_circle_center_to_end(1)),2*pi);
    is_counterClockwise = fcn_geometry_arcDirectionFrom3Points(this_island_points(min_index,:), this_island_points(best_SNR_index,:), this_island_points(max_index,:),-1);
    best_arc_flag_is_counterclockwise  = (1==is_counterClockwise);
    best_arc_flag_this_is_a_circle     = 0;

    best_fit_arc = [...
        best_arc_center,...
        best_arc_radius,...
        best_arc_start_angle_in_radians,...
        best_arc_end_angle_in_radians,...
        best_arc_flag_this_is_a_circle,...
        best_arc_flag_is_counterclockwise,...
        ];

    % Plot the arcs?
    if flag_do_debug
        figure(debug_fig_num)
        subplot(1,3,1);
        fcn_geometry_plotGeometry('arc',best_fit_arc);
    end

    % Save results
    best_fit_arcs = [best_fit_arcs; best_fit_arc]; %#ok<AGROW>

end

% Plot final results?
if ~isempty(fig_num)
    figure(fig_num);
    clf;

    hold on;
    grid on;
    axis equal
    xlabel('X [m]');
    ylabel('Y [m]');

    % Plot the results
    for ith_arc = 1:length(best_fit_arcs(:,1))
        fcn_geometry_plotGeometry('arc',best_fit_arcs(ith_arc,:));
    end

    % Plot the input points
    plot(this_island_points(:,1),this_island_points(:,2),'b.','MarkerSize',5);

    % Set the axis range
    temp = axis;
    axis_range_x = temp(2)-temp(1);
    axis_range_y = temp(4)-temp(3);
    percent_larger = 0.3;
    axis([temp(1)-percent_larger*axis_range_x, temp(2)+percent_larger*axis_range_x,  temp(3)-percent_larger*axis_range_y, temp(4)+percent_larger*axis_range_y]);

end
end % Ends fcn_INTERNAL_curvatureArcsFromSNR


%% fcn_INTERNAL_curvatureModelSort
function [fitSequence_fitTypes, fitSequence_parameters, model_SNRs, sorted_model_fit_ID_at_each_index, arc_matrix, segment_matrix] = fcn_INTERNAL_curvatureModelSort(revised_good_arcs, all_segments, best_fit_SNRs, model_fit_ID_at_each_index, fig_num)

% How many points do we have?
Npoints = length(model_fit_ID_at_each_index(:,1));

% Initialize output cell arrays and matricies
fitSequence_parameters = {};
fitSequence_fitTypes   = {};
arc_matrix = [];
segment_matrix = [];
NumFits = 0;

% Initialize the variable storing the model ordering
model_ordering = [];

% Find the first model
index_of_first_model = find(model_fit_ID_at_each_index~=0,1);

current_fit_type = 0;
for ith_point = index_of_first_model:Npoints
    this_fit = model_fit_ID_at_each_index(ith_point);
    if 0~=this_fit
        if this_fit~=current_fit_type
            NumFits = NumFits+1;
            current_fit_type = this_fit;
            model_ordering = [model_ordering; this_fit]; %#ok<AGROW>
            if current_fit_type>0
                fitSequence_fitTypes{NumFits}   = 'arc'; %#ok<AGROW>
                fitSequence_parameters{NumFits} = revised_good_arcs(current_fit_type,:); %#ok<AGROW>
                arc_matrix = [arc_matrix; revised_good_arcs(current_fit_type,:)]; %#ok<AGROW>
            else
                fitSequence_fitTypes{NumFits}   = 'segment'; %#ok<AGROW>
                fitSequence_parameters{NumFits} = all_segments(current_fit_type*-1,:); %#ok<AGROW>
                segment_matrix = [segment_matrix; all_segments(current_fit_type*-1,:)]; %#ok<AGROW>
            end
        end
    end
end

% Fix the model ordering and assign SNRs
sorted_model_fit_ID_at_each_index = 0*model_fit_ID_at_each_index;
for ith_model = 1:length(model_ordering)
    original_model_number = model_ordering(ith_model);
    indicies_to_update = model_fit_ID_at_each_index==original_model_number;
    sorted_model_fit_ID_at_each_index(indicies_to_update,1) = ith_model;
end

% Check if SNRs need to be updated
if ~isempty(best_fit_SNRs)
    model_SNRs = 0*best_fit_SNRs;
    for ith_model = 1:length(model_ordering)
        original_model_number = model_ordering(ith_model);
        model_SNRs(ith_model,1) = best_fit_SNRs(original_model_number,1);
    end
else
    model_SNRs = [];
end

% Plot the domain fits
Nmodels = length(fitSequence_fitTypes);
for ith_model = 1:Nmodels
    color = [1 0 0]*(1 - (ith_model-1)/Nmodels);
    fcn_geometry_plotGeometry(fitSequence_fitTypes{ith_model},fitSequence_parameters{ith_model},[],color,fig_num);
end
   
end % Ends fcn_INTERNAL_curvatureModelSort


%% fcn_INTERNAL_alignArcsBySNRC2
function arc_matrix_C2 = fcn_INTERNAL_alignArcsBySNRandC2(arc_matrix, model_SNRs, best_fit_number_at_each_index, fig_num)
% Given fitSequence_fitTypes that are all arcs, and the arc's
% fitSequence_parameters with a SNR for each arc fit, goes through the arcs
% in order of SNR from highest to lowest, making sure that each arc is C2
% continuous with all its adjacent arcs

Npoints = length(best_fit_number_at_each_index);

% Make sure the entries are same length
NumArcs = length(arc_matrix(:,1));

% NumArcs = length(fitSequence_fitTypes);
% assert(length(fitSequence_parameters)==NumArcs);
assert(length(model_SNRs)==NumArcs);

% % Make sure all entries are 'arc' types
% assert(all(strcmp(fitSequence_fitTypes,'arc')));
% 
% % Save all the arcs into a matrix, as it's easier to work with
% arc_matrix = zeros(NumArcs,length(fitSequence_parameters{1}(1,:)));
% for ith_arc_ranking = 1:NumArcs
%     arc_matrix(ith_arc_ranking,:) = fitSequence_parameters{ith_arc_ranking}(1,:);
% end

% Sort the arcs by SNR
[~,sort_order] = sort(model_SNRs,'descend');
[~,inverse_sort_order] = sort(sort_order);

transverse_tolerance = 2;

% Initialize the output arc parameter array
arcs_with_C2_continuity = arc_matrix;

fitSequence_fitTypes{NumArcs} = [];
fitSequence_parameters{NumArcs} = [];

for ith_arc = 1:NumArcs
    fitSequence_fitTypes{ith_arc} = 'arc';
    fitSequence_parameters{ith_arc} = arc_matrix(ith_arc,:);
end
fcn_geometry_plotFitSequences(fitSequence_fitTypes, fitSequence_parameters,[],[0.2 0.2 0.2],(fig_num));

% Proceed through the arcs from highest SNR to lowest
for ith_arc_ranking = 1:NumArcs
    this_arc = sort_order(ith_arc_ranking);
    fcn_geometry_plotGeometry('arc',arc_matrix(this_arc,:),[],[1 0 0],fig_num);

    if ith_arc_ranking>1

        if this_arc>1
            adjacent_before_model = this_arc-1;
        else
            adjacent_before_model = [];
        end
        if this_arc<NumArcs
            adjacent_after_model = this_arc+1;
        else
            adjacent_after_model = [];
        end

        % Only worry about rankings that are lower than this one
        adjacent_models = [adjacent_before_model; adjacent_after_model];
        rankings_to_check = inverse_sort_order(adjacent_models);
        rankings_to_check = rankings_to_check(rankings_to_check<ith_arc_ranking);
                
        % Take the minimum of all the adjacent rankings. If all adjacent
        % rankings are Nan, this returns NaN
        [priority_ranking, index_match] = min(rankings_to_check);

        % Get the lowest ranked (e.g. best) arc that is adjacent, and use
        % this as the reference
        flag_is_feasible = 1;
        if ~isempty(priority_ranking) && ~isnan(priority_ranking)
            adjacent_model_to_prioritize = adjacent_models(index_match);
            [flag_is_feasible, ~, closest_feasible_arc2_parameters] = ...
                fcn_geometry_isC2FeasibleArcToArc(arc_matrix(adjacent_model_to_prioritize,:), arc_matrix(this_arc,:), (transverse_tolerance), (0.0001), (3456));
            arcs_with_C2_continuity(this_arc,:) = closest_feasible_arc2_parameters;
        end

        % Was it feasible?
        if 1~=flag_is_feasible
            % Method that would fix this: perform C1 continuity, then shift
            % C1 line segment to get C2 continuity. 
            error('Unable to merge arcs');
            % % Perform C1 continuity 
            % [revised_arc1_parameters, revised_arc2_parameters, revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters]  = ...
            %     % fcn_geometry_alignArcArc(arc1_parameters, arc2_parameters, (threshold), (continuity_level),  (fig_num))

        end

        % Check the other side - if there's a curve there also, need to
        % make sure it still fits. If it doesn't fit, do not move it - this
        % would break the code.
        if flag_is_feasible && length(rankings_to_check)==2

            if index_match==1
                secondary_match = 2;
            else
                secondary_match = 1;
            end
            adjacent_model_to_prioritize = adjacent_models(secondary_match);

            [flag_is_feasible, ~, closest_feasible_arc2_parameters] = ...
                fcn_geometry_isC2FeasibleArcToArc(arc_matrix(adjacent_model_to_prioritize,:), closest_feasible_arc2_parameters, transverse_tolerance, (0.0001), (3456));

        end

        if 1~=flag_is_feasible
            error('Unable to fit with C2 curvature')
            % Method to fix this: remove this arc from the list and attempt
            % to connect the two adjacent to each other with C2 and then C1
        end
    end
end

arc_matrix_C2 = arcs_with_C2_continuity; 

fitSequence_fitTypes_C2 = fitSequence_fitTypes;
fitSequence_parameters_C2 = fitSequence_parameters; 
for ith_arc_ranking = 1:NumArcs
    fitSequence_parameters_C2{ith_arc_ranking} = arcs_with_C2_continuity(ith_arc_ranking,:);
end

% Plot the revised results
figure(fig_num);
clf;

hold on;
grid on;
axis equal
xlabel('X [m]');
ylabel('Y [m]');

% Plot the input sequence:
format_string = sprintf(' ''-'',''Color'',[0.6 0.6 0.6],''LineWidth'',4 ');
fcn_geometry_plotFitSequences(fitSequence_fitTypes, fitSequence_parameters,[],format_string,(fig_num));

% Plot the output sequence
format_string = sprintf(' ''-'',''Color'',[0 0 1],''LineWidth'',2 ');
fcn_geometry_plotFitSequences(fitSequence_fitTypes_C2, fitSequence_parameters_C2,[],format_string,(fig_num));

end % Ends fcn_INTERNAL_alignArcsBySNRC2
