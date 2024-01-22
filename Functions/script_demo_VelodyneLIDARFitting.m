% script_demo_VelodyneLIDARFitting

clear all; %#ok<CLALL>
close all;


%% Load the data
load('PointCloud_Separated_Data.mat')
load('PointCloud_Separated_Data_Raw.mat');

% Load the RingNames and N_rings
currentScan = ptCloud_pts_layers_separated_cell{1};
currentScanRaw = ptCloud_pts_raw_rings_cell_save{1};
RingNames = fieldnames(currentScan);
N_rings = length(RingNames);

%% For plotting

% Get the color ordering?
try
    color_ordering = orderedcolors('gem12');
catch
    color_ordering = colororder;
end

N_colors = length(color_ordering(:,1));

%% Look at data
flag_do_plots = 0;
good_limits = [0.5    2.5    1.0    3   -1    1];

% time_iteration = 44; % 2, 4, 5, 6, 7*, 10, 23*, 32, 35*, 44, 45*

for time_iteration = 1:length(ptCloud_pts_layers_separated_cell)
    % BAD ones: 8 is bad, 9 is mediocre, 11 is hard, 12 should work but doesnt
    % 24 is a good test case

    %  The radius of the sphere target is 0.151 meter, and all the points are in meters

    % Columns 1:3 are XYZ in LIDAR coordinates
    % 4 is intensity
    % 5 in Ring number
    % 6 is the time difference (milliseconds?)

    % The following if statement is only to plot the data, and to update
    % "good_limits" for plotting
    if flag_do_plots
        figure(4575);
        clf;
        hold on;
        grid on;
        axis equal;



        currentScan = ptCloud_pts_layers_separated_cell{time_iteration};
        currentScanRaw = ptCloud_pts_raw_rings_cell_save{time_iteration};


        % Plot just the separated scans, so we can get the axis
        for ith_ring = 1:N_rings
            ringName = RingNames{ith_ring};
            RingData = currentScan.(ringName);
            plot3(RingData(:,1), RingData(:,2), RingData(:,3),'.','MarkerSize',10);
            if ~isempty(RingData)
                text(RingData(1,1), RingData(1,2), RingData(1,3), sprintf('%.0d',ith_ring));
            end
        end
        view(3)

        good_limits = axis;

        % Now plot the entire scan
        for ith_ring = 1:N_rings
            ringName = RingNames{ith_ring};
            RingData = currentScan.(ringName);
            RingDataRaw = currentScanRaw.(ringName);
            plot3(RingDataRaw(:,1), RingDataRaw(:,2), RingDataRaw(:,3),'k.','MarkerSize',10);
            plot3(RingData(:,1), RingData(:,2), RingData(:,3),'.','MarkerSize',30);
            if ~isempty(RingData)
                text(RingData(1,1), RingData(1,2), RingData(1,3), sprintf('%.0d',ith_ring));
            end
        end
        view(3)
        axis(good_limits);

    end


    % Pull out scans
    currentScan = ptCloud_pts_layers_separated_cell{time_iteration};
    currentScanRaw = ptCloud_pts_raw_rings_cell_save{time_iteration};

    if flag_do_plots
        fig_num = 10;
        figure(fig_num);
        clf;
        hold on;
        grid on;
        axis equal;

        fig_num = 20;
        figure(fig_num);
        clf;
        hold on;
        grid on;
        axis equal;
    end

    % Set constraints on Hough fit
    transverse_tolerance = 0.02;
    station_tolerance = 3;
    points_required_for_agreement = 10;
    flag_force_circle_fit = 1;
    nominalRadius = 0.151; % meters
    expected_radii_range = [0.09 0.2];
    flag_use_permutations = [];

    % Try Hough fitting, then regression fitting, to pull out all the centers
    xyz_centers = [];
    true_xyz_centers = [];
    rings_to_keep = [];

    N_domains = N_rings;
    Hough_domains_to_keep{N_domains} = struct;
    Regression_domains_to_keep{N_domains} = struct;

    for ith_ring = 1:N_rings  % 2:6
        ringName = RingNames{ith_ring};
        RingData = currentScan.(ringName);
        dataToFit = RingData(:,1:2);

        if ~isempty(dataToFit)
            % Perform Hough fit
            fig_num = 10;
            inputPoints = dataToFit;
            if length(inputPoints(:,1))>points_required_for_agreement
                Hough_domains  = ...
                    fcn_geometry_fitHoughCircle(inputPoints, transverse_tolerance, ...
                    (station_tolerance), (points_required_for_agreement), (flag_force_circle_fit), (expected_radii_range), (flag_use_permutations), -1);

                if flag_do_plots
                    fcn_geometry_plotFitDomains(Hough_domains, fig_num);
                end

                % Save best Hough results
                Hough_domains_to_keep{ith_ring} = Hough_domains{1};

                % Check the regression fit
                fig_num = 20;
                regression_domains = fcn_geometry_HoughRegression(Hough_domains, -1);

                % Save best fit results
                Regression_domains_to_keep{ith_ring} = Hough_domains{1};


                if ~isnan(regression_domains{1}.best_fit_parameters)
                    mean_z = mean(RingData(:,3)); % Replace this line with calculated z-height
                    circleCenter = regression_domains{1}.best_fit_parameters(1:2);
                    circleRadius = regression_domains{1}.best_fit_parameters(3);
                    if flag_do_plots
                        fcn_geometry_plotFitDomains(regression_domains, fig_num);
                        fprintf(1,'Best-fit radius for ring %.0d: %.3f meters\n', ith_ring, circleRadius);
                    end

                    if 1==1
                        xyz_centers = [xyz_centers; circleCenter mean_z]; %#ok<AGROW>
                        rings_to_keep = [rings_to_keep; ith_ring]; %#ok<AGROW>
                    else

                        % Calculate height error - the error in z-height to the
                        % center of the sphere, based on how small the measured
                        % radius is to the true radius
                        height_error = 0;
                        if circleRadius < nominalRadius
                            height_error = (nominalRadius.^2 - circleRadius.^2).^0.5;

                            % There are 2 possible votes for XYZ of center
                            xyz_centers = [xyz_centers; circleCenter mean_z+height_error]; %#ok<AGROW>
                            rings_to_keep = [rings_to_keep; ith_ring]; %#ok<AGROW>
                            xyz_centers = [xyz_centers; circleCenter mean_z-height_error]; %#ok<AGROW>
                            rings_to_keep = [rings_to_keep; ith_ring]; %#ok<AGROW>

                        else
                            xyz_centers = [xyz_centers; circleCenter mean_z]; %#ok<AGROW>
                            rings_to_keep = [rings_to_keep; ith_ring]; %#ok<AGROW>
                            true_xyz_centers = [true_xyz_centers; circleCenter mean_z]; %#ok<AGROW>
                        end
                    end

                end
            end

        end
    end

    if flag_do_plots
        figure(4575);
        plot3(xyz_centers(:,1), xyz_centers(:,2), xyz_centers(:,3),'.-','MarkerSize',30, 'Linewidth',3,'Color',[0 0 1]);

        figure(2343);
        plot(xyz_centers(:,1), xyz_centers(:,2),'.','MarkerSize',30, 'Linewidth',3,'Color',[0 0 1]);
    end



    % Keep only the circles close to each other
    distance_center_to_center_threshold = 0.1;
    flags_count_of_ring_is_close = zeros(length(rings_to_keep(:,1)),1);
    for ith_ring = 1:length(rings_to_keep)
        distances = sum((xyz_centers(ith_ring,1:2)-xyz_centers(:,1:2)).^2,2).^0.5;
        good_distances = (distances<distance_center_to_center_threshold);
        count_of_close_distances = sum(good_distances)-1; % Subtract the self-count
        flags_count_of_ring_is_close(ith_ring) = count_of_close_distances;
    end

    max_agreement = max(flags_count_of_ring_is_close);
    indicies_of_rings_that_agree = find(flags_count_of_ring_is_close==max_agreement);
    rings_in_agreement = rings_to_keep(indicies_of_rings_that_agree);

    % Find the center of the circle
    % mean_circle_center = 0.5*mean(xyz_centers(indicies_of_rings_that_agree,:)) + 0.5*mean(true_xyz_centers);
    mean_circle_center = mean(xyz_centers(indicies_of_rings_that_agree,:),1);

    if flag_do_plots
        figure(4575);
        plot3(mean_circle_center(:,1), mean_circle_center(:,2), mean_circle_center(:,3),'r.','MarkerSize',50, 'Linewidth',3);
    end

    % Grab all the good domains
    N_good_domains = length(rings_in_agreement);
    good_Hough_domains_to_keep{N_good_domains} = struct;
    good_Regression_domains_to_keep{N_good_domains} = struct;
    for ith_ring = 1:length(rings_in_agreement)
        ring_to_keep = rings_in_agreement(ith_ring);
        good_Hough_domains_to_keep{ith_ring} = Hough_domains_to_keep{ring_to_keep};
        good_Regression_domains_to_keep{ith_ring} = Regression_domains_to_keep{ring_to_keep};
    end

    if flag_do_plots
        fcn_geometry_plotFitDomains(good_Hough_domains_to_keep, 2345);
    end

    % Save results
    results_to_plot(time_iteration).mean_circle_center = mean_circle_center; %#ok<SAGROW>
    results_to_plot(time_iteration).rings_in_agreement = rings_in_agreement;  %#ok<SAGROW>
    results_to_plot(time_iteration).good_Regression_domains_to_keep = good_Regression_domains_to_keep;  %#ok<SAGROW>

end % Ends loop through time
save("Data\LIDAR_debugging.mat","results_to_plot");

%% Load and plot the data
load("Data\LIDAR_debugging.mat","results_to_plot");
load('PointCloud_Separated_Data_Raw.mat'); % Contains ptCloud_pts_raw_rings_cell_save

% Load the RingNames and N_rings
% currentScan = ptCloud_pts_layers_separated_cell{1};
currentScanRaw = ptCloud_pts_raw_rings_cell_save{1};
RingNames = fieldnames(currentScanRaw);
N_rings = length(RingNames);

% For plotting
% Get axis limits
good_limits = [0.5    2.5    1.0    3   -1    1];


% Get the color ordering
try
    color_ordering = orderedcolors('gem12');
catch
    color_ordering = colororder;
end
N_colors = length(color_ordering(:,1));


for time_iteration = 1:95 % length(ptCloud_pts_raw_rings_cell_save)
    currentScanRaw = ptCloud_pts_raw_rings_cell_save{time_iteration};

    mean_circle_center = results_to_plot(time_iteration).mean_circle_center;
    rings_in_agreement = results_to_plot(time_iteration).rings_in_agreement;  
    good_Regression_domains_to_keep = results_to_plot(time_iteration).good_Regression_domains_to_keep;  

    fcn_INTERNAL_plotResults(time_iteration, mean_circle_center, RingNames, currentScanRaw, rings_in_agreement, good_Regression_domains_to_keep, good_limits, color_ordering);
    pause(0.01);
    
end % Ends loop through time

%% Find and plot all the points in each ring, in 3D
function fcn_INTERNAL_plotResults(time_iteration, mean_circle_center, RingNames, currentScanRaw, rings_in_agreement, good_Regression_domains_to_keep, good_limits, color_ordering)
figure(47464);
clf;
hold on;
grid on;
axis equal

title(sprintf('Time: %.0d',time_iteration));

plot3(mean_circle_center(:,1), mean_circle_center(:,2), mean_circle_center(:,3),'r.','MarkerSize',50, 'Linewidth',3);


N_colors = length(color_ordering(:,1));

% Start by plotting the entire scan
N_rings = length(RingNames);
for ith_ring = 1:N_rings
    ringName = RingNames{ith_ring};
    RingDataRaw = currentScanRaw.(ringName);
    plot3(RingDataRaw(:,1), RingDataRaw(:,2), RingDataRaw(:,3),'k.','MarkerSize',10);
end
view(30,30)
axis(good_limits);



% Initialize arrays
all_close_RingData = {};
all_fitted_RingData = {};
only_fitted_sphere_points = [];
only_sphere_points = [];

if 1==0
    % Pull out and plot all the fits
    for ith_ring = 1:length(rings_in_agreement)
        current_ring = rings_in_agreement(ith_ring);
        ringName = RingNames{current_ring};
        RingDataRaw = currentScanRaw.(ringName);
        pointData3d = RingDataRaw(:,1:3);

        % Find the fitted data
        % pointData2d = good_Hough_domains_to_keep{ith_ring}.points_in_domain;
        pointData2d = good_Regression_domains_to_keep{ith_ring}.points_in_domain;
        fitted_RingData = [];
        for ith_point = 1:length(pointData2d(:,1))
            % Find where the points match to 5 decimal places
            N_decimals = 4;
            point2d = floor(pointData2d(ith_point,1:2)*(10^N_decimals));
            points3d = floor(pointData3d(:,1:2)*(10^N_decimals));
            [flag_match_find, index_of_match] = ismember(point2d, points3d, 'rows');
            if ~flag_match_find
                error('no match found');
            end

            fitted_RingData = [fitted_RingData; RingDataRaw(index_of_match,:)]; %#ok<AGROW>
        end

        % Save the results
        all_fitted_RingData{ith_ring} = fitted_RingData; %#ok<AGROW>


        % Plot results
        current_color = color_ordering(mod(ith_ring,N_colors)+1,:);
        plot3(all_fitted_RingData{ith_ring}(:,1), all_fitted_RingData{ith_ring}(:,2), all_fitted_RingData{ith_ring}(:,3),'.','MarkerSize',30,'Color',current_color);
    end
end

% Pull out and plot all the close points in any ring.
% Set a distance threshold
threshold_point_distance = 0.2;
for ith_ring = 1:N_rings
    ringName = RingNames{ith_ring};
    RingDataRaw = currentScanRaw.(ringName);
    pointData3d = RingDataRaw(:,1:3);

    % Find the close data
    distance_to_mean_center = sum((pointData3d(:,1:3)-mean_circle_center(1,1:3)).^2,2).^0.5;
    close_data = RingDataRaw(distance_to_mean_center<threshold_point_distance,:);

    % Save the results
    all_close_RingData{ith_ring} = close_data; %#ok<AGROW>
    only_fitted_sphere_points = [only_fitted_sphere_points; close_data]; %#ok<AGROW>

    % Plot results
    plot3(all_close_RingData{ith_ring}(:,1), all_close_RingData{ith_ring}(:,2), all_close_RingData{ith_ring}(:,3),'.','MarkerSize',15,'Color',[1 0 0]);
end

plot3(only_fitted_sphere_points(:,1), only_fitted_sphere_points(:,2), only_fitted_sphere_points(:,3),'.','MarkerSize',10,'Color',[1  1 0]);
end



%% Grab the equations for the planes of each scan

% fig_num = 2221;
%
%
% currentScan = ptCloud_pts_layers_separated_cell{1};
% RingNames = fieldnames(currentScan);
%
% % Loop through the rings1
% for ith_ring = 1:N_rings
%     ringName = RingNames{ith_ring};
%
%     for time_iteration = 1:10 %length(ptCloud_pts_layers_separated_cell)
%         currentScan = ptCloud_pts_layers_separated_cell{time_iteration};
%         currentScanRaw = ptCloud_pts_raw_rings_cell_save{time_iteration};
%
%         RingData = currentScan.(ringName);
%         RingDataRaw = currentScanRaw.(ringName);
%
%         if ~isempty(RingData)
%             % Plot the data
%             figure(fig_num);
%             clf;
%             hold on;
%             grid on;
%             axis equal;
%             view(3);
%
%
%             plot3(RingDataRaw(:,1), RingDataRaw(:,2), RingDataRaw(:,3),'k.','MarkerSize',10);
%             plot3(RingData(:,1), RingData(:,2), RingData(:,3),'.','MarkerSize',30);
%             text(RingData(1,1), RingData(1,2), RingData(1,3), sprintf('%.0d',ith_ring));
%             axis(good_limits);
%
%             % Find the plane fit
%             [fitted_parameters, standard_deviation_in_z, z_fit] = fcn_geometry_fitPlaneLinearRegression(RingData(:,1:3),123);
%             axis(good_limits);
%
%             figure(37372)
%             clf;
%             hold on;
%             grid on;
%             axis equal;
%
%             plot(RingData(:,3),z_fit,'.');
%
%             % Plot the best-fit line
%             temp = axis;
%             min_value = min(temp);
%             max_value = max(temp);
%             plot([min_value max_value],[min_value max_value],'k-');
%             axis(temp);
%
%             % Plot x versus z
%             xy_fig = 38383;
%             figure(xy_fig);
%             clf;
%             hold on;
%             grid on;
%             axis equal;
%
%             % Do the rawdata fitting
%             plot(RingData(:,1),RingData(:,3),'k.','MarkerSize',20)
%             plot(RingData(:,1),z_fit,'g.','MarkerSize',10)
%             temp_axis = axis;
%
%             [rawfitted_parameters, rawstandard_deviation_in_z, z_fit] = fcn_geometry_fitPlaneLinearRegression(RingDataRaw(:,1:3),456);
%             axis(good_limits);
%
%             figure(37373)
%             clf;
%             hold on;
%             grid on;
%             axis equal;
%
%             plot(RingDataRaw(:,3),z_fit,'.');
%
%
%             % Plot the best-fit line
%             temp = axis;
%             min_value = min(temp);
%             max_value = max(temp);
%             plot([min_value max_value],[min_value max_value],'k-');
%             axis(temp);
%
%
%             figure(xy_fig);
%             plot(RingDataRaw(:,1),RingDataRaw(:,3),'r.','MarkerSize',20)
%             plot(RingDataRaw(:,1),z_fit,'g.','MarkerSize',10)
%             axis(temp_axis)
%
%
%             fprintf(1,'Ring: %.0d, time: %.0d\n',ith_ring,time_iteration);
%             params = fitted_parameters;
%             fprintf(1,'\tFit of RANSAC  : %.5f \t %.5f \t %.5f \t with STD: %.5f meters \n',params(1),params(2),params(3),standard_deviation_in_z);
%             params = rawfitted_parameters;
%             fprintf(1,'\tFit of raw     : %.5f \t %.5f \t %.5f \t with STD: %.5f meters \n',params(1),params(2),params(3),standard_deviation_in_z);
%
%             % Check means
%             mean_x = mean(RingData(:,1));
%             mean_y = mean(RingData(:,2));
%             plane_z = [mean_x mean_y 1]*fitted_parameters; % Solve for z vertices data
%             plane_z_raw = [mean_x mean_y 1]*rawfitted_parameters; % Solve for z vertices data
%             fprintf(1,'\tZ-heights (meters) RANSAC: %.5f \n',plane_z);
%             fprintf(1,'\tZ-heights (meters) RAW   : %.5f \n',plane_z_raw);
%
%             pause(0.01);
%
%         end
%     end
%
% end
