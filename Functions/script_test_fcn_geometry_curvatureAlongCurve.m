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



%% Calculate curvatures to find any islands in geometric information
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
if Npoints>100
    data_width = ceil(minimum_island_separation/average_spacing);
    if data_width<=2
        data_width = [];
        warning('Not enough data for curvature calculations. Using all data.');
    end
    % 20 points is usually plenty
    data_width = max(data_width,20);
end

fig_num = 1234;


[curvatures, arc_centers, index_ranges, point_curvature_minimum, curvature_SNRs] = fcn_geometry_curvatureAlongCurve(points_to_fit, (data_width), (fig_num));

% Assign islands to locations where the SNR is less than 1, e.g. it's more
% likely that the data is a line than an arc. We make it 3 here to give it
% a bit of wiggle-room (some are right on edge).
is_island = curvature_SNRs>3;

%% Fix loops
% In the test track, and in any data that forms a loop, it is common that
% the start/end of data will be on a curve. The results will be an island
% at start and at end, but with straightaways somehwere before and after
% the end. 
%
% If this is the case, they can be fixed by shifting the data to start at
% the first "lake" part after the island.

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

%% Find islands
% TO DO: recode this as a for loop, not while loop. The number of islands
% is equal to the number of transitions from 0 to 1 in the is_island vector

N_islands = 0;
island_ranges = {};
filled_islands = is_island;

while any(filled_islands)

    island_starts = find(filled_islands,1);


    % Make everything an island up to where the island starts. This makes
    % finding where the island ends very trivial
    filled_islands(1:island_starts) = 1;
    island_ends = find(filled_islands==0,1);

    % If nothing found, then everything remaining was an island
    if isempty(island_ends)
        island_ends = Npoints;
    end

    % Set everything up to end of island a lake
    filled_islands(1:island_ends) = 0;

    % Update count and information for the island
    N_islands = N_islands+1;
    island_ranges{N_islands} = (island_starts:island_ends)'; %#ok<SAGROW>


end

% Plot the results
fit_fig_num = 67543;
figure(fit_fig_num);
clf;

hold on;
grid on;
axis equal
xlabel('X [m]');
ylabel('Y [m]');

plot(shifted_points(:,1),shifted_points(:,2),'k.','MarkerSize',20);
for ith_island = 1:N_islands
    this_island_range = island_ranges{ith_island};
    plot(shifted_points(this_island_range,1),shifted_points(this_island_range,2),'.','MarkerSize',10);
end

%% Calculate the curvatures of each island
fig_num = 2346;

for ith_island = 1:1 %N_islands
    this_island_range = island_ranges{ith_island};
    this_island_points = shifted_points(this_island_range,:);

    data_width = []; % Default is to use all possible points
    [curvatures, arc_centers, index_ranges, point_curvature_minimum, curvature_SNRs] = fcn_geometry_curvatureAlongCurve(this_island_points, (data_width), (fig_num+ith_island));

end


%% Pull out the fits by ranking them via SNRs
N_fits = 0;
best_fit_arcs = [];
best_fit_SNRs = [];
best_fit_ranges = [];
remaining_curvature_SNR = curvature_SNRs;

remaining_curvature_SNR(remaining_curvature_SNR<20) = nan;

best_fit_number_at_each_index = nan*curvature_SNRs;

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
        min_index = best_SNR_index-index_ranges(best_SNR_index);
        max_index = best_SNR_index+index_ranges(best_SNR_index);
        this_index_range = (min_index:max_index)';
        plot(this_island_points(min_index:max_index,1),this_island_points(min_index:max_index,2),'m.','MarkerSize',10)

        % Plot the max SNR point
        plot(this_island_points(best_SNR_index,1),this_island_points(best_SNR_index,2),'g.','MarkerSize',30)

        axis(temp_axis);

        title('Input points');

        %%%%%%%%%%%%%%%%%%%%%%
        subplot(1,3,2);
        % cla;

        semilogy(this_island_range,curvatures,'k-');
        hold on;
        semilogy(this_island_range,point_curvature_minimum,'-','Color',[0.6 0.6 0.6]);

        % Plot the part covered by this fit
        plot(this_island_range(this_index_range),curvatures(best_SNR_index)*ones(length(this_index_range)),'k-','Markersize',20,'LineWidth',3);

        grid on;
        xlabel('index [count]');
        ylabel('curvature [1/m]');
        title('Curvatures')

        %%%%%%%%%%%%%%%%%%%%%%
        subplot(1,3,3);
        % cla;
        grid on;
        hold on;

        plot(this_island_range, remaining_curvature_SNR,'k-');

        % Plot the part covered by this fit
        plot(this_island_range(this_index_range), remaining_curvature_SNR(this_index_range,:),'k-','LineWidth',3);
        plot(this_island_range(best_SNR_index), remaining_curvature_SNR(best_SNR_index),'g.','Markersize',30);

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
    this_index_range = (min_index:max_index)';

    remaining_curvature_SNR(this_index_range) = nan;
    is_zeros = find(best_fit_number_at_each_index(this_index_range)==0);
    unfilled_index_range = this_index_range(is_zeros);
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


%% Plot raw results

SNR_threshold = 30;
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
