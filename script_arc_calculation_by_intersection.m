%% script_test_fcn_geometry_fitSequentialArcsViaRegions
% Exercises the function: fcn_geometry_fitSequentialArcns

% 2024_04_14 - S. Brennan
% -- wrote the code


rng(1)

%% Use fillArcSequence to create some test data
fig_num = 1;
figure(fig_num);
clf;

rng(1); % Fix the random number, for debugging

% arc_pattern has [1/R and L] for each segment as a row
% arc_pattern = [...
%     1/20, 15; 
%     0 20;
%     -1/5 10; 
%     0 10;
%     1/15 40; 
%     0 15
%     -1/10 20];

arc_pattern = [...
    1/10, 15; 
    0 20];

M = 10;
sigma = 0.02;

[test_points, ~, ~, trueArcStartIndicies, trueNamedCurveTypes, trueParameters] = fcn_geometry_fillArcSequenceTestPoints(arc_pattern, M, sigma, -1);

% Add noise?
if 1==0
    % Corrupt the results
    probability_of_corruption = 1;
    magnitude_of_corruption = 0.03;

    test_points = fcn_geometry_corruptPointsWithOutliers(test_points,...
        (probability_of_corruption), (magnitude_of_corruption), (-1));
end

% Add outliers?
if 1==0
    % Corrupt the results
    probability_of_corruption = 0.1;
    magnitude_of_corruption = 1;

    test_points = fcn_geometry_corruptPointsWithOutliers(test_points,...
        (probability_of_corruption), (magnitude_of_corruption), (-1));
end

fig_num = 1;
figure(fig_num);
clf;

fitting_tolerance = 0.1; % Units are meters
maxRadius = 300; %the maximum possible radius

%%%% 
% Step 0, find the first point more than the tolerance away from the
% starting point
far_enough_index = fcn_INTERNAL_findNextTestPointIndex(test_points,fitting_tolerance,1);
far_enough_point = test_points(far_enough_index,:);
indicies_points_missed = (2:(far_enough_index-1))';

% Now find the point after the last missed point that is also more than the
% tolerance away, so we can "turn around" later.
cleared_missed_points_index = fcn_INTERNAL_findNextTestPointIndex(test_points,fitting_tolerance,points_missed(end));

if 1~=flag_found_point_outside_search
    error('Unable to find any point more than the search tolerance away from the start point. Cannot continue.');
end

%%%%
debug_fig_num = [];
if ~isempty(debug_fig_num)
    figure(debug_fig_num); clf;
end
% Step 1, find bounds for starting points
[old_cumulative_positive_region, old_cumulative_negative_region] = fcn_INTERNAL_findTangentVectorsForPoints(starting_point,far_enough_point,fitting_tolerance, maxRadius, debug_fig_num);

%%%%
% Step 2, find the regions that clear the 1st point up to cleared the
% missed point
flag_found_final_region = 0;
for ith_point = far_enough_index+1:cleared_missed_points_index
    [positive_region, negative_region] = fcn_INTERNAL_findTangentVectorsForPoints(test_points(1,:),test_points(ith_point,:),fitting_tolerance, maxRadius, debug_fig_num);

    if old_cumulative_positive_region.NumRegions>0
        cumulative_positive_region = intersect(old_cumulative_positive_region,positive_region);
    end
    if old_cumulative_negative_region.NumRegions>0
        cumulative_negative_region = intersect(old_cumulative_negative_region,negative_region);
    end
    if cumulative_negative_region.NumRegions==0 && cumulative_positive_region.NumRegions==0
        flag_found_final_region = 1;
        break;
    end


    old_cumulative_positive_region = cumulative_positive_region;
    old_cumulative_negative_region = cumulative_negative_region;

    figure(fig_num);
    clf;
    hold on

    plot(trueParameters{1}(1,1),trueParameters{1}(1,2),'r.','Markersize',40);
    plot(test_points(1:ith_point,1),test_points(1:ith_point,2),'k.','MarkerSize',10);
    plot(cumulative_positive_region);
    plot(cumulative_negative_region);

    pause(0.01);
end

%%%%
% Step 3, find the regions for the missed points
if 0==flag_found_final_region
    for ith_point = points_missed(1):points_missed(end)
        [positive_region, negative_region] = fcn_INTERNAL_findTangentVectorsForPoints(test_points(ith_point,:),test_points(cleared_missed_points_index,:),fitting_tolerance, maxRadius, debug_fig_num);

        if old_cumulative_positive_region.NumRegions>0
            cumulative_positive_region = intersect(old_cumulative_positive_region,positive_region);
        end
        if old_cumulative_negative_region.NumRegions>0
            cumulative_negative_region = intersect(old_cumulative_negative_region,negative_region);
        end
        if cumulative_negative_region.NumRegions==0 && cumulative_positive_region.NumRegions==0
            flag_found_final_region = 1;
            break;
        end


        old_cumulative_positive_region = cumulative_positive_region;
        old_cumulative_negative_region = cumulative_negative_region;

        figure(fig_num);
        clf;
        hold on

        plot(trueParameters{1}(1,1),trueParameters{1}(1,2),'r.','Markersize',40);
        plot(test_points(1:ith_point,1),test_points(1:ith_point,2),'k.','MarkerSize',10);
        plot(cumulative_positive_region);
        plot(cumulative_negative_region);

        pause(0.01);

    end
end

%%%%
% Step 4, find the regions after cleared missed point

if 0==flag_found_final_region
    for ith_point = cleared_missed_points_index(end)+1:NtestPoints
        [positive_region, negative_region] = fcn_INTERNAL_findTangentVectorsForPoints(test_points(1,:),test_points(ith_point,:),fitting_tolerance, maxRadius, debug_fig_num);

        if old_cumulative_positive_region.NumRegions~=0
            cumulative_positive_region = intersect(old_cumulative_positive_region,positive_region);
            if isequal(old_cumulative_positive_region,cumulative_positive_region)
                %flag_found_final_region = 1;
                numfails = numfails+1;
            else
                numfails = 0;
            end
                
        end
        if old_cumulative_negative_region.NumRegions~=0
            cumulative_negative_region = intersect(old_cumulative_negative_region,negative_region);
            if isequal(old_cumulative_negative_region,cumulative_negative_region)
                %flag_found_final_region = 1;
                numfails = numfails+1;
            else
                numfails = 0;
            end
        end

        if numfails >=3
            % break;
        end

        uncertain_positive = polybuffer(cumulative_positive_region,-2*fitting_tolerance);
        uncertain_negative = polybuffer(cumulative_negative_region,-2*fitting_tolerance);


        if ith_point>174
            figure(74747);
            clf;
            hold on;
            plot(trueParameters{1}(1,1),trueParameters{1}(1,2),'r.','Markersize',40);
            plot(old_cumulative_positive_region);
            plot(positive_region);
            plot(cumulative_positive_region);
            plot(uncertain_positive);

        end


        if (uncertain_positive.NumRegions==0) && (uncertain_negative.NumRegions==0)
            flag_found_final_region = 1;
        end

        if flag_found_final_region
            break;
        end


        old_cumulative_positive_region = cumulative_positive_region;
        old_cumulative_negative_region = cumulative_negative_region;

        figure(fig_num);
        clf;
        hold on

        plot(trueParameters{1}(1,1),trueParameters{1}(1,2),'r.','Markersize',40);
        plot(test_points(1:ith_point,1),test_points(1:ith_point,2),'k.','MarkerSize',10);
        plot(cumulative_positive_region);
        plot(cumulative_negative_region);

        pause(0.01);


    end
end

if old_cumulative_positive_region.NumRegions>0
    final_region = old_cumulative_positive_region;
end
if old_cumulative_negative_region.NumRegions>0
    final_region = old_cumulative_negative_region;
end

[x,y] = centroid(final_region)

% Compare to the other fitting method
% [-0.0943 20.2800 20.2881 4.7170 5.5595 0 1]

flag_fit_backwards = 0;
[fitSequence_points_forward, fitSequence_shapes_forward, fitSequence_endIndicies_forward, fitSequence_parameters_forward, fitSequence_bestFitType_forward] = ...
    fcn_geometry_fitSequentialArcs(test_points, fitting_tolerance, flag_fit_backwards, [], 2);

fitSequence_parameters_forward{1}(1:2)

% [distance, location, path_segment, t, u] = ...
%         fcn_Path_findProjectionHitOntoPath(path,...
%         sensor_vector_start,sensor_vector_end,...
%         (4),(fig_num));



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

%%
function far_enough_index = fcn_INTERNAL_findNextTestPointIndex(test_points,fitting_tolerance,starting_index)
% Finds the index of the next point that is more than 2 times the fitting
% tolerance away from the starting_point, given by a starting index.
starting_point = test_points(starting_index,:);
next_test_point_index = starting_index+1;
flag_found_point_outside_search = 0;
NtestPoints = length(test_points(:,1));

far_enough_index = nan;
while (next_test_point_index<=NtestPoints) && (0==flag_found_point_outside_search)
    far_enough_point = test_points(next_test_point_index,:);
    distance_this_point_to_start = abs(sum((starting_point - far_enough_point).^2,2).^0.5);
    if distance_this_point_to_start > 2*fitting_tolerance
        flag_found_point_outside_search = 1;
        far_enough_index = next_test_point_index;
    end
    next_test_point_index = next_test_point_index + 1;
end

end

%% fcn_INTERNAL_findTangentVectorsForPoints
function [positive_region, negative_region] = fcn_INTERNAL_findTangentVectorsForPoints(point_1,point_2,fitting_tolerance, maxRadius, fig_num)
centers_start = point_1;
centers_end   = point_2;
radii_start   = fitting_tolerance;
radii_end     = fitting_tolerance;
flag_inside_or_outside = -1;
voting_points_start = [];
voting_points_end   = [];

[points_tangent_start, points_tangent_end] = ...
    fcn_geometry_findTangentPointsTwoCircles(...
    centers_start,...
    centers_end,...
    radii_start,...
    radii_end,...
    flag_inside_or_outside,...
    voting_points_start,voting_points_end,...
    fig_num);

% Find and plot the centerpoint
centerPoint = (point_1+point_2)/2;

% Figure out which vector is the bottom one relative to the start circle.
% The bottom vector is the one where the cross-product direction is
% positive when crossing the radius vector on the starting circle toward
% the center point.
start_circle_radius_vector = points_tangent_start(1,:) - point_1;
start_circle_to_center_vector = centerPoint - point_1;
cross_product_result = cross([start_circle_radius_vector 0],[start_circle_to_center_vector 0]);
if cross_product_result(3)>0
    index_low_vector_on_start  = 1;
    index_high_vector_on_start = 2;
else
    index_low_vector_on_start  = 2;
    index_high_vector_on_start = 1;
end

Low1High2_vector = points_tangent_end(index_low_vector_on_start,:) - points_tangent_start(index_low_vector_on_start,:);
High1Low2_vector = points_tangent_end(index_high_vector_on_start,:) - points_tangent_start(index_high_vector_on_start,:);
unit_High1Low2_vector = fcn_geometry_calcUnitVector(High1Low2_vector);
unit_Low1High2_vector = fcn_geometry_calcUnitVector(Low1High2_vector);

unit_ortho_High1Low2_vector = unit_High1Low2_vector*[0 1; -1 0];
unit_ortho_Low1High2_vector = unit_Low1High2_vector*[0 1; -1 0];

center_vector = (unit_ortho_High1Low2_vector+unit_ortho_Low1High2_vector)/2;
unit_center_vector = fcn_geometry_calcUnitVector(center_vector);

points_positive_region = ones(4,1)*centerPoint + maxRadius*[0 0;  unit_ortho_High1Low2_vector; unit_center_vector;  unit_ortho_Low1High2_vector];
points_negative_region = ones(4,1)*centerPoint - maxRadius*[0 0;  unit_ortho_High1Low2_vector; unit_center_vector;  unit_ortho_Low1High2_vector];

positive_region = polyshape(points_positive_region(:,1),points_positive_region(:,2));
negative_region = polyshape(points_negative_region(:,1),points_negative_region(:,2));

if ~isempty(fig_num)
    figure(fig_num)
    plot(centerPoint(:,1),centerPoint(:,2),'k.','MarkerSize',10);

    quiver(centerPoint(1,1),centerPoint(1,2),unit_ortho_High1Low2_vector(1,1),unit_ortho_High1Low2_vector(1,2),0,'Color','b');
    quiver(centerPoint(1,1),centerPoint(1,2),unit_ortho_Low1High2_vector(1,1),unit_ortho_High1Low2_vector(1,2),0,'Color','g');
    quiver(centerPoint(1,1),centerPoint(1,2),unit_center_vector(1,1),unit_center_vector(1,2),0,'Color','r');

    plot(positive_region);
    plot(negative_region);
end
end % Ends fcn_INTERNAL_findTangentVectorsForPoints

