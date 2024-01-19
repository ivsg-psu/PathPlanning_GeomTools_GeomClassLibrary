% script_demo_VelodyneLIDARFitting

clear all;
close all;


%% Load the data
load('PointCloud_Separated_Data.mat')

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
RingNames = fieldnames(currentScan);

for ith_ring = 1:length(RingNames)
    ringName = RingNames{ith_ring};
    RingData = currentScan.(ringName);
    plot3(RingData(:,1), RingData(:,2), RingData(:,3),'.','MarkerSize',30);
    if ~isempty(RingData)
        text(RingData(1,1), RingData(1,2), RingData(1,3), sprintf('%.0d',ith_ring));
    end
end
view(3)

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


