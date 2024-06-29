%% script_test_fcn_geometry_curvatureAlongCurve
% Exercises the function: fcn_geometry_curvatureAlongCurve

% 2024_06_27 - S. Brennan
% -- wrote the code

close all;

%% Basic example - normal case in arc

%
fig_num = 1;
figure(fig_num);
clf;

rng(1); % Fix the random number, for debugging


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
    keep_every = 100; % 20 works OK
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

%%
% Code starts here

Npoints = length(points_to_fit(:,1));
data_width = []; % 50;

% Initialize storage arrays
index_curvatures = nan(Npoints,1);
curvatures = nan(Npoints,1);
arc_centers = nan(Npoints,2);
index_ranges = nan(Npoints,1);
point_curvature_minimums = nan(Npoints,1);

% For each point, calculate it's best-SNR curvature, and save result
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

%% Pull out the fits by ranking them via SNRs
N_fits = 0;
best_fit_arcs = [];
best_fit_SNRs = [];
best_fit_ranges = [];
remaining_curvature_SNR = curvature_SNR;

flag_do_debug = 1;
debug_fig_num = 38383;
figure(debug_fig_num);clf;

flag_first_time = 1;

while ~all(isnan(remaining_curvature_SNR))

    % Find the best remaining
    [~,best_SNR_index] = max(remaining_curvature_SNR);
    
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

        plot(best_SNR_index,point_curvature_minimums(best_SNR_index),'g.','Markersize',20);

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
    remaining_curvature_SNR(min_index:max_index) = nan;

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


%% Plot results
% URHERE
% best_SNR_index = 410;

%%%%%%%%%%%%%%%%%%%%%%
figure(fig_num)
subplot(1,3,1);
cla;

hold on;
grid on;
axis equal
xlabel('X [m]');
ylabel('Y [m]');

% Plot the input points
plot(points_to_fit(:,1),points_to_fit(:,2),'b.','MarkerSize',20);

% Make axis slightly larger?
temp = axis;
axis_range_x = temp(2)-temp(1);
axis_range_y = temp(4)-temp(3);
percent_larger = 0.3;
axis([temp(1)-percent_larger*axis_range_x, temp(2)+percent_larger*axis_range_x,  temp(3)-percent_larger*axis_range_y, temp(4)+percent_larger*axis_range_y]);

temp_axis = axis;

% Plot the circle fit at the point
fcn_geometry_plotCircle(arc_centers(best_SNR_index,:), 1/curvatures(best_SNR_index),'r-',fig_num);


% Plot the index range
min_index = best_SNR_index-index_ranges(best_SNR_index);
max_index = best_SNR_index+index_ranges(best_SNR_index);
plot(points_to_fit(min_index:max_index,1),points_to_fit(min_index:max_index,2),'m.','MarkerSize',10)

% Plot the max SNR point
plot(points_to_fit(best_SNR_index,1),points_to_fit(best_SNR_index,2),'g.','MarkerSize',30)

axis(temp_axis);

title('Input points');

%%%%%%%%%%%%%%%%%%%%%%
subplot(1,3,2);
cla;

semilogy(index_curvatures,curvatures,'k-');
hold on;
semilogy(index_curvatures,point_curvature_minimums,'-','Color',[0.6 0.6 0.6]);

plot(best_SNR_index,point_curvature_minimums(best_SNR_index),'g.','Markersize',20);

grid on;
xlabel('index [count]');
ylabel('curvature [1/m]');
title('Curvatures')

%%%%%%%%%%%%%%%%%%%%%%
subplot(1,3,3);
cla;
grid on;
hold on;

plot(index_curvatures, curvature_SNR,'k-');
plot(index_curvatures(best_SNR_index), curvature_SNR(best_SNR_index),'g.','Markersize',30);

xlabel('index [count]');
ylabel('SNR [unitless]');
title('Curvature SNR')

% Check sizes
assert(isequal(size(point_curvature),[1 1]));
assert(isequal(size(point_circle_center),[1 2]));
assert(isequal(size(index_range),[1 1]));
assert(isequal(size(point_curvature_minimum),[1 1]));

%% Test track test

