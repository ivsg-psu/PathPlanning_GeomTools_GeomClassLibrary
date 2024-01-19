% script_demo_VelodyneLIDARFitting

clear all;
close all;


%% Load the data
load('PointCloud_Separated_Data.mat')
load('PointCloud_Separated_Data_Raw.mat');

%% Look at data
time_iteration = 10; % 2, 4, 5, 6, 7*, 10
% BAD ones: 8 is bad, 9 is mediocre

%  The radius of the sphere target is 0.151 meter, and all the points are in meters

% Columns 1:3 are XYZ in LIDAR coordinates
% 4 is intensity
% 5 in Ring number
% 6 is the time difference (milliseconds?)

figure(4575);
clf;
hold on;
grid on;
axis equal;


currentScan = ptCloud_pts_layers_separated_cell{time_iteration};
currentScanRaw = ptCloud_pts_raw_rings_cell_save{time_iteration};
RingNames = fieldnames(currentScan);

% Plot just the separated scans, so we can get the axis
for ith_ring = 1:length(RingNames)
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
for ith_ring = 1:length(RingNames)
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

%% Grab the equations for the planes of each scan
fig_num = 2221;


currentScan = ptCloud_pts_layers_separated_cell{1};
RingNames = fieldnames(currentScan);

% Loop through the rings1
for ith_ring = 1:length(RingNames)
    ringName = RingNames{ith_ring};

    for time_iteration = 1:10 %length(ptCloud_pts_layers_separated_cell)
        currentScan = ptCloud_pts_layers_separated_cell{time_iteration};
        currentScanRaw = ptCloud_pts_raw_rings_cell_save{time_iteration};

        RingData = currentScan.(ringName);
        RingDataRaw = currentScanRaw.(ringName);

        if ~isempty(RingData)
            % Plot the data
            figure(fig_num);
            clf;
            hold on;
            grid on;
            axis equal;
            view(3);


            plot3(RingDataRaw(:,1), RingDataRaw(:,2), RingDataRaw(:,3),'k.','MarkerSize',10);
            plot3(RingData(:,1), RingData(:,2), RingData(:,3),'.','MarkerSize',30);
            text(RingData(1,1), RingData(1,2), RingData(1,3), sprintf('%.0d',ith_ring));
            axis(good_limits);

            % Find the plane fit
            [fitted_parameters, standard_deviation_in_z, z_fit] = fcn_geometry_fitPlaneLinearRegression(RingData(:,1:3),123);
            axis(good_limits);

            figure(37372)
            clf;
            hold on;
            grid on;
            axis equal;

            plot(RingData(:,3),z_fit,'.');

            % Plot the best-fit line
            temp = axis;
            min_value = min(temp);
            max_value = max(temp);
            plot([min_value max_value],[min_value max_value],'k-');
            axis(temp);

            % Plot x versus z
            xy_fig = 38383;
            figure(xy_fig);
            clf;
            hold on;
            grid on;
            axis equal;
            
            % Do the rawdata fitting
            plot(RingData(:,1),RingData(:,3),'k.','MarkerSize',20)
            plot(RingData(:,1),z_fit,'g.','MarkerSize',10)
            temp_axis = axis;

            [rawfitted_parameters, rawstandard_deviation_in_z, z_fit] = fcn_geometry_fitPlaneLinearRegression(RingDataRaw(:,1:3),456);
            axis(good_limits);

            figure(37373)
            clf;
            hold on;
            grid on;
            axis equal;

            plot(RingDataRaw(:,3),z_fit,'.');


            % Plot the best-fit line
            temp = axis;
            min_value = min(temp);
            max_value = max(temp);
            plot([min_value max_value],[min_value max_value],'k-');
            axis(temp);


            figure(xy_fig);
            plot(RingDataRaw(:,1),RingDataRaw(:,3),'r.','MarkerSize',20)
            plot(RingDataRaw(:,1),z_fit,'g.','MarkerSize',10)
            axis(temp_axis)


            fprintf(1,'Ring: %.0d, time: %.0d\n',ith_ring,time_iteration);
            params = fitted_parameters;
            fprintf(1,'\tFit of RANSAC  : %.5f \t %.5f \t %.5f \t with STD: %.5f meters \n',params(1),params(2),params(3),standard_deviation_in_z);
            params = rawfitted_parameters;
            fprintf(1,'\tFit of raw     : %.5f \t %.5f \t %.5f \t with STD: %.5f meters \n',params(1),params(2),params(3),standard_deviation_in_z);
            
            % Check means
            mean_x = mean(RingData(:,1));
            mean_y = mean(RingData(:,2));
            plane_z = [mean_x mean_y 1]*fitted_parameters; % Solve for z vertices data
            plane_z_raw = [mean_x mean_y 1]*rawfitted_parameters; % Solve for z vertices data
            fprintf(1,'\tZ-heights (meters) RANSAC: %.5f \n',plane_z);
            fprintf(1,'\tZ-heights (meters) RAW   : %.5f \n',plane_z_raw);
            
            pause(0.01);

        end
    end

end

%% Pull out scans
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

% Set parameters on Hough fit
transverse_tolerance = 0.02;
station_tolerance = 3;
points_required_for_agreement = 10;
flag_force_circle_fit = 1;
expected_radii_range = [0.09 0.2];
flag_use_permutations = [];

% Try fitting
xyz_centers = [];
for ith_ring = 2:6 % 1:length(RingNames)
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
                (station_tolerance), (points_required_for_agreement), (flag_force_circle_fit), (expected_radii_range), (flag_use_permutations), []);
            fcn_geometry_plotFitDomains(Hough_domains, fig_num);

            % Check the regression fit
            fig_num = 20;
            regression_domains = fcn_geometry_HoughRegression(Hough_domains, -1);

            if ~isnan(regression_domains{1}.best_fit_parameters)
                mean_z = mean(RingData(:,3));
                circleCenter = regression_domains{1}.best_fit_parameters(1:2);
                circleRadius = regression_domains{1}.best_fit_parameters(3);
                fcn_geometry_plotFitDomains(regression_domains, fig_num);
                fprintf(1,'Best-fit radius for ring %.0d: %.3f meters\n', ith_ring, circleRadius);
                xyz_centers = [xyz_centers; circleCenter mean_z]; %#ok<AGROW>
            end
        end

    end
end

figure(4575);
plot3(xyz_centers(:,1), xyz_centers(:,2), xyz_centers(:,3),'.-','MarkerSize',30, 'Linewidth',3,'Color',[0 0 1]);


