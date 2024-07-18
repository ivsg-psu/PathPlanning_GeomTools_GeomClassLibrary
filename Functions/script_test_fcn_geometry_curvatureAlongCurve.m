%% script_test_fcn_geometry_curvatureAlongCurve
% Exercises the function: fcn_geometry_curvatureAlongCurve

% 2024_06_27 - S. Brennan
% -- wrote the code

close all;

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
    keep_every = 1; % 100 works OK
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
    %     -1/5 10;
    %     0 10];

    M = 1; % How many points per meter
    sigma = 0.02; % The standard deviation in the points relative to the perfect function fit, in meters

    [points_to_fit, ~, ~, trueArcStartIndicies, trueNamedCurveTypes, trueParameters] = fcn_geometry_fillArcSequenceTestPoints(arc_pattern, M, sigma, -1);

end

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

spatial_differences = diff(points_to_fit(:,1:2));
spatial_distances   = real(sum(spatial_differences.^2,2).^0.5);
average_spacing     = mean(spatial_distances);

% Check if the data is a loop (this changes a few steps that follow)
distance_XY_from_start = real(sum((points_to_fit(:,1:2) - points_to_fit(1,1:2)).^2,2).^0.5);
max_distance_XY_from_start = max(distance_XY_from_start);
distance_start_to_end  = real(sum((points_to_fit(end,1:2) - points_to_fit(1,1:2)).^2,2).^0.5);
flag_is_a_loop = 0;
if distance_start_to_end<(5*average_spacing)
    flag_is_a_loop = 1;
end

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




%% Calculate any islands in geometric information
% Islands are where there are interconnected arcs that have no line
% segments within. Islands are useful for analysis because calculations
% done within an island are not affected by calculations in other islands,
% and so the data can be sub-grouped by island and thereby saving huge
% amounts of computation.
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
minimum_island_separation = 1; % Units are meters
data_width = ceil(minimum_island_separation/average_spacing); 
if data_width<=2
    error('Not enough data for curvature calculations. Quitting.');
end


% Initialize storage arrays
index_curvatures = nan(Npoints,1); % Records which index was tested
curvatures = nan(Npoints,1);  % Records the best-fit curvature at each point
arc_centers = nan(Npoints,2); % Records the arc centers of each point
index_ranges = nan(Npoints,1); % Records how "wide" the curvature calculation was able to expand
point_curvature_minimums = nan(Npoints,1); % Records the largest radius (min curvature) that could possibly fit the curvature (best case)

% For each point, calculate curvature and SNR, and save result
for ith_point = 1:Npoints
    fprintf(1,'Testing point %.0d of %.0d\n',ith_point,Npoints);
    [point_curvature, point_circle_center, index_range, point_curvature_minimum] = fcn_geometry_curvatureAtPoint(points_to_fit, ith_point, data_width, (-1));
    index_curvatures(ith_point,1) = ith_point;
    curvatures(ith_point,1) = point_curvature;
    arc_centers(ith_point,:) = point_circle_center;
    index_ranges(ith_point,1) = index_range;
    point_curvature_minimums(ith_point,1) = point_curvature_minimum;
end

% Fix NaN values at start/end of data. These are areas where the curvature
% is undefined since curves cannot be fit at the end of the data
% (regression requires minimum 3 points). To fix, we have the end points
% inheret the values from their closest neighbors.

N_to_fix = find(~isnan(curvatures),1)-1;

% Fix NaN values at start
index_to_inheret = N_to_fix+1;
range_to_fix = 1:N_to_fix;
curvatures(range_to_fix,1) = curvatures(index_to_inheret,1);
arc_centers(range_to_fix,1) = arc_centers(index_to_inheret,1);
arc_centers(range_to_fix,2) = arc_centers(index_to_inheret,2);
index_ranges(range_to_fix,1) = index_ranges(index_to_inheret,1);
point_curvature_minimums(range_to_fix,1) = point_curvature_minimums(index_to_inheret,1);

% Fix NaN values at end
index_to_inheret = Npoints - N_to_fix;
range_to_fix = (Npoints - N_to_fix + 1):Npoints;
curvatures(range_to_fix,1) = curvatures(index_to_inheret,1);
arc_centers(range_to_fix,1) = arc_centers(index_to_inheret,1);
arc_centers(range_to_fix,2) = arc_centers(index_to_inheret,2);
index_ranges(range_to_fix,1) = index_ranges(index_to_inheret,1);
point_curvature_minimums(range_to_fix,1) = point_curvature_minimums(index_to_inheret,1);



% Plot the results
%%%%%%%%%%%%%%%%%%%%%%
island_curvature_fig_num = 444;
figure(island_curvature_fig_num)
clf;

subplot(1,2,1);

semilogy(index_curvatures,curvatures,'k-');
hold on;
semilogy(index_curvatures,point_curvature_minimums,'-','Color',[0.6 0.6 0.6]);

grid on;
xlabel('index [count]');
ylabel('curvature [1/m]');
title('Curvatures')

%%%%%%%%%%%%%%%%%%%%%%
subplot(1,2,2);

semilogy(index_curvatures, curvature_SNR,'k-');
hold on;

xlabel('index [count]');
ylabel('SNR [unitless]');
title('Curvature SNR')


% Assign islands to locations where the SNR is less than 1, e.g. it's more
% likely that the data is a line than an arc
is_island = curvature;




%% Calculate the curvatures at each point
% Code starts here

data_width = []; % The number of data points to consider to right and left of the test point. Default is to use them all

% Initialize storage arrays
index_curvatures = nan(Npoints,1);
curvatures = nan(Npoints,1);
arc_centers = nan(Npoints,2);
index_ranges = nan(Npoints,1);
point_curvature_minimums = nan(Npoints,1);

% For each point, calculate curvature and SNR, and save result
for ith_point = 1:Npoints
    fprintf(1,'Testing point %.0d of %.0d\n',ith_point,Npoints);
    [point_curvature, point_circle_center, index_range, point_curvature_minimum] = fcn_geometry_curvatureAtPoint(points_to_fit, ith_point, data_width, (-1));
    index_curvatures(ith_point,1) = ith_point;
    curvatures(ith_point,1) = point_curvature;
    arc_centers(ith_point,:) = point_circle_center;
    index_ranges(ith_point,1) = index_range;
    point_curvature_minimums(ith_point,1) = point_curvature_minimum;
end

% Fix line data, which can have bad curvatures
point_curvature_minimums(curvatures<0.001) = 1;

curvature_SNR = curvatures./point_curvature_minimums;

%% Plot the results
fit_fig_num = 2343;
figure(fit_fig_num);
clf;

hold on;
grid on;
axis equal
xlabel('X [m]');
ylabel('Y [m]');

plot(points_to_fit(:,1),points_to_fit(:,2),'k.','MarkerSize',20);
bad_points = find(curvature_SNR<20);
plot(points_to_fit(bad_points,1),points_to_fit(bad_points,2),'r.','MarkerSize',10);


%% Pull out the fits by ranking them via SNRs
N_fits = 0;
best_fit_arcs = [];
best_fit_SNRs = [];
best_fit_ranges = [];
remaining_curvature_SNR = curvature_SNR;

remaining_curvature_SNR(remaining_curvature_SNR<20) = nan;

best_fit_number_at_each_index = nan*curvature_SNR;

flag_do_debug = 1;
debug_fig_num = 38383;
figure(debug_fig_num);clf;

flag_first_time = 1;

NumFits = 0;

while ~all(isnan(remaining_curvature_SNR))

    % Find the best remaining
    [~,best_SNR_index] = max(remaining_curvature_SNR);
    NumFits = NumFits+1;
    
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
            plot(points_to_fit(:,1),points_to_fit(:,2),'b.','MarkerSize',20);


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
        min_index = best_SNR_index-index_ranges(best_SNR_index);
        max_index = best_SNR_index+index_ranges(best_SNR_index);
        index_range = (min_index:max_index)';
        plot(points_to_fit(min_index:max_index,1),points_to_fit(min_index:max_index,2),'m.','MarkerSize',10)

        % Plot the max SNR point
        plot(points_to_fit(best_SNR_index,1),points_to_fit(best_SNR_index,2),'g.','MarkerSize',30)

        axis(temp_axis);

        title('Input points');

        %%%%%%%%%%%%%%%%%%%%%%
        subplot(1,3,2);
        % cla;

        semilogy(index_curvatures,curvatures,'k-');
        hold on;
        semilogy(index_curvatures,point_curvature_minimums,'-','Color',[0.6 0.6 0.6]);

        % Plot the part covered by this fit
        plot(index_range,curvatures(best_SNR_index)*ones(length(index_range)),'k-','Markersize',20,'LineWidth',3);

        grid on;
        xlabel('index [count]');
        ylabel('curvature [1/m]');
        title('Curvatures')

        %%%%%%%%%%%%%%%%%%%%%%
        subplot(1,3,3);
        % cla;
        grid on;
        hold on;

        plot(index_curvatures, remaining_curvature_SNR,'k-');

        % Plot the part covered by this fit
        plot(index_range, remaining_curvature_SNR(index_range,:),'k-','LineWidth',3);
        plot(index_curvatures(best_SNR_index), remaining_curvature_SNR(best_SNR_index),'g.','Markersize',30);

        xlabel('index [count]');
        ylabel('SNR [unitless]');
        title('Curvature SNR')
    end
    
    % Save results
    best_fit_SNRs    = [best_fit_SNRs; remaining_curvature_SNR(best_SNR_index)]; %#ok<AGROW>
    best_fit_ranges = [best_fit_ranges; index_ranges(best_SNR_index)]; %#ok<AGROW>
    

    % Block out the indicies of this fit
    min_index = best_SNR_index-index_ranges(best_SNR_index);
    max_index = best_SNR_index+index_ranges(best_SNR_index);
    index_range = (min_index:max_index)';

    remaining_curvature_SNR(index_range) = nan;
    is_zeros = find(best_fit_number_at_each_index(index_range)==0);
    unfilled_index_range = index_range(is_zeros);
    best_fit_number_at_each_index(unfilled_index_range) = NumFits;

    % Save fit parameters
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
    vector_from_circle_center_to_start = points_to_fit(min_index,:)-best_arc_center;
    vector_from_circle_center_to_end   = points_to_fit(max_index,:)-best_arc_center;
    best_arc_start_angle_in_radians    = mod(atan2(vector_from_circle_center_to_start(2),vector_from_circle_center_to_start(1)),2*pi);
    best_arc_end_angle_in_radians      = mod(atan2(vector_from_circle_center_to_end(2),vector_from_circle_center_to_end(1)),2*pi);
    is_counterClockwise = fcn_geometry_arcDirectionFrom3Points(points_to_fit(min_index,:), points_to_fit(best_SNR_index,:), points_to_fit(max_index,:),-1);
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


%% Plot raw results

SNR_threshold = 20;
good_fits = find(best_fit_SNRs>SNR_threshold);
good_arcs = best_fit_arcs(good_fits,:);
NumFitsGood = length(good_fits);


%% Ensure C2 continuity
transverse_tolerance = 2;
revised_good_arcs = good_arcs;

fit_fig_num = 2343;
figure(fit_fig_num);
clf;

hold on;
grid on;
axis equal
xlabel('X [m]');
ylabel('Y [m]');

for ith_fit = 1:length(good_arcs(:,1))
    fcn_geometry_plotGeometry('arc',good_arcs(ith_fit,:));
    if ith_fit>1

        % Find any prior fits that are adjacent
        first_index = find(best_fit_number_at_each_index==ith_fit,1,'first');
        last_index  = find(best_fit_number_at_each_index==ith_fit,1,'last');
        
        if first_index>1 && best_fit_number_at_each_index(first_index-1)<ith_fit
            adjacent_to_first = best_fit_number_at_each_index(first_index-1);
        else
            adjacent_to_first = nan;
        end

        if last_index<Npoints && best_fit_number_at_each_index(last_index+1)<ith_fit
            adjacent_to_last = best_fit_number_at_each_index(last_index+1);
        else
            adjacent_to_last = nan;
        end
        data_to_check = [adjacent_to_first adjacent_to_last];

        [priority_match, index_match] = min(data_to_check);

        if ~isnan(priority_match)
            [flag_is_feasible, feasibility_distance, closest_feasible_arc2_parameters] = ...
                fcn_geometry_isC2FeasibleArcToArc(best_fit_arcs(priority_match,:), best_fit_arcs(ith_fit,:), (transverse_tolerance/2), (0.0001), (3456));
            revised_good_arcs(ith_fit,:) = closest_feasible_arc2_parameters;
        end

        % Check the other side - if there's a curve there also, need to
        % make sure it still fits
        if flag_is_feasible && ~any(isnan(data_to_check))

            if index_match==1
                secondary_match =data_to_check(2);
            else
                secondary_match = data_to_check(1);
            end
            [flag_is_feasible, feasibility_distance, closest_feasible_arc2_parameters] = ...
                fcn_geometry_isC2FeasibleArcToArc(best_fit_arcs(secondary_match,:), closest_feasible_arc2_parameters, 0, (0.0001), (3456));

        end            

        if 1~=flag_is_feasible
            error('Unable to fit with C2 curvature')
        end
    end
end

%% Plot the revised results
final_fig_num = 6456;
figure(final_fig_num);
clf;

hold on;
grid on;
axis equal
xlabel('X [m]');
ylabel('Y [m]');

for ith_fit = 1:length(good_arcs(:,1))
    fcn_geometry_plotGeometry('arc',revised_good_arcs(ith_fit,:));
end

%% Fill the unfilled areas
indicies_unfilled = best_fit_number_at_each_index>NumFitsGood;

all_segments = [];

num_segments = 0;
while any(indicies_unfilled)
    num_segments = num_segments+1;

    start_index = find(indicies_unfilled==1,1);
    remainder   = indicies_unfilled;
    remainder(1:start_index) = 1; % Fill in any zeros beforehand
    end_index   = find(remainder==0,1); % Find first zero afterwards

    % Shut off the indicies searched thus far, in prep for next round
    indicies_unfilled(1:end_index) = 0;

    % Check for errors
    if 1==start_index
        error('Start index should never be 1');
    end
    if isempty(end_index)
        error('Empty end indicies should never occur')
    end
    
    % Grab the arc before the open area
    fit_before = best_fit_number_at_each_index(start_index-1);
    fit_after  = best_fit_number_at_each_index(end_index);
    arc_before = revised_good_arcs(fit_before,:);
    arc_after  = revised_good_arcs(fit_after,:);

    continuity_level = 1;
    [revised_arc1_parameters, revised_arc2_parameters, revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters]  = ...
        fcn_geometry_alignArcArc(arc_before, arc_after, (transverse_tolerance), (continuity_level),  (57894));

    % Save arc revisions
    revised_good_arcs(fit_before,:) = revised_arc1_parameters;
    revised_good_arcs(fit_after,:) = revised_arc2_parameters;
    
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
    best_fit_number_at_each_index(start_index:(end_index-1),1) = -1*num_segments;

end

%% Plot the revised results
arc_line_fig_num = 23434;
figure(arc_line_fig_num);
clf;

hold on;
grid on;
axis equal
xlabel('X [m]');
ylabel('Y [m]');

for ith_fit = 1:length(good_arcs(:,1))
    fcn_geometry_plotGeometry('arc',revised_good_arcs(ith_fit,:));
end


for ith_fit = 1:length(all_segments(:,1))
    fcn_geometry_plotGeometry('segment',all_segments(ith_fit,:));
end

%% Turn this into a sequence of fits
fitSequence_parameters = {};
fitSequence_fitTypes   = {};

current_fit_type = best_fit_number_at_each_index(3);
NumFits = 1;

if current_fit_type>0
    fitSequence_fitTypes{NumFits}   = 'arc';
    fitSequence_parameters{NumFits} = revised_good_arcs(current_fit_type,:);
else
    fitSequence_fitTypes{NumFits}   = 'segment';
    fitSequence_parameters{NumFits} = all_segments(current_fit_type*-1,:);
end


for ith_point = 3:Npoints
    this_fit = best_fit_number_at_each_index(ith_point);
    if ~isnan(this_fit)
        if this_fit~=current_fit_type
            NumFits = NumFits+1;
            current_fit_type = this_fit;
            if current_fit_type>0
                fitSequence_fitTypes{NumFits}   = 'arc';
                fitSequence_parameters{NumFits} = revised_good_arcs(current_fit_type,:);
            else
                fitSequence_fitTypes{NumFits}   = 'segment';
                fitSequence_parameters{NumFits} = all_segments(current_fit_type*-1,:);
            end
        end
    end
end

% Plot the domain fits
fig_num = 1111;
fcn_geometry_plotFitSequences(fitSequence_fitTypes, fitSequence_parameters,[],[],(fig_num));

%% Align results
fig_num = 76767;
continuity_level = 2;
revised_fitSequence_parameters = ...
    fcn_geometry_alignGeometriesInSequence(fitSequence_fitTypes, fitSequence_parameters, transverse_tolerance, (continuity_level), (fig_num));
