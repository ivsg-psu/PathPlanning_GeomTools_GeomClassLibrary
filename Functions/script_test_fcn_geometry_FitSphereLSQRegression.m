% script_test_fcn_geometry_FitSphereLSQRegression.m
% This is a script to exercise the function: fcn_geometry_FitSphereLSQRegression.m
% This function was written on 2023_11_22 by Xinyu Cao, xfc5113@psu.edu

% Revision history:
% 2023_11_22
% -- first write of the code


%% Set up the workspace
close all
clc

%% Basic example
fig_num = 1;

% Fill in the true values


% Call the fitting function
[C_sphere,R_sphere,E_total] = fcn_geometry_FitSphereLSQRegression(XYZ_array);


%% ADVANCED EXAMPLE
% Use the GPSandLidar_Data.mat in the Data folder to start, the loaded data
% is a struct containing two fields
% rawdata_selected.GPS_SparkFun_Temp_ENU_selected and
% rawdata_selected.Lidar_pointcloud_cell_selected
load('GPSandLidar_Data.mat','rawdata_selected');

%% Grab the GPS and Lidar Data
GPS_SparkFun_Temp_ENU = rawdata_selected.GPS_SparkFun_Temp_ENU_selected;
Lidar_pointcloud_cell = rawdata_selected.Lidar_pointcloud_cell_selected;

%% Create Fake PointCloud Data
N_fake_points = 20;
thetavec = linspace(0,pi,N_fake_points).';
phivec = linspace(0,2*pi,2*N_fake_points).';
[theta_fake, phi_fake] = meshgrid(thetavec, phivec);


% Simulate a sphere using the same radius of the sphere target
R_set = 0.151;
r_vec = R_set*ones(size(phi_fake));
N_data = length(GPS_SparkFun_Temp_ENU);
gps_select_range = 1:N_data;
noise_magnitude = 0.01; % the magnitude of the random noise added to the data
fake_PointCloud_cell{N_data} = [];
for n_gps_idxs = gps_select_range
    center = GPS_SparkFun_Temp_ENU(n_gps_idxs,:);
    X_fake = center(1) + r_vec.*sin(theta_fake).*cos(phi_fake);
    Y_fake = center(2) + r_vec.*sin(theta_fake).*sin(phi_fake);
    Z_fake = center(3) + r_vec.*cos(theta_fake);
    size_fake_data = size(X_fake);
    X_fake_reshaped = reshape(X_fake,[size_fake_data(1)*size_fake_data(2),1]);
    Y_fake_reshaped = reshape(Y_fake,[size_fake_data(1)*size_fake_data(2),1]);
    Z_fake_reshaped = reshape(Z_fake,[size_fake_data(1)*size_fake_data(2),1]);
    
    % Generate random noise to X, Y and Z
    X_noise = noise_magnitude*randn(size(X_fake_reshaped));
    Y_noise = noise_magnitude*randn(size(Y_fake_reshaped));
    Z_noise = noise_magnitude*randn(size(Z_fake_reshaped));
    
    % Add the noise to X, Y and Z
    X_fake_with_noise = X_fake_reshaped+X_noise;
    Y_fake_with_noise = Y_fake_reshaped+Y_noise;
    Z_fake_with_noise = Z_fake_reshaped+Z_noise;
    XYZ_array_with_noise = [X_fake_with_noise, Y_fake_with_noise, Z_fake_with_noise];
    Intensity_fake = round(255*rand(size(X_fake_reshaped)));
    fake_PointCloud_cell{n_gps_idxs,1} = [XYZ_array_with_noise,Intensity_fake];
end

%% View fake point cloud data
Lidar_origin = mean(GPS_SparkFun_Temp_ENU,1);
x_range = [Lidar_origin(1)-1.5, Lidar_origin(1)+1.5];
y_range = [Lidar_origin(2)-1.5, Lidar_origin(2)+1.5];
z_range = [Lidar_origin(3)-1.5, Lidar_origin(3)+1.5];
pause_time = 0.01; % pause 0.01 second
LidarViewer = pcplayer(x_range,y_range,z_range,'ColorSource','Intensity');
for i_scan = gps_select_range
    current_scan = fake_PointCloud_cell{i_scan,1};
    XYZ_array = current_scan(:,1:3);
    Intensity = current_scan(:,4);
    ptCloud = pointCloud(XYZ_array, 'Intensity', Intensity);
    if ~isnan(ptCloud.Location)
        view(LidarViewer,ptCloud)
        LidarViewer.Axes.Title.String = sprintf("Frame number: %d", i_scan);
        pause(pause_time)
    end
end

%% Fit Sphere Target with the fake pointcloud data
clear C_array R_array E_ave_array pointCloud_space_cell

for idx_fit = 1:N_data
    current_scan = fake_PointCloud_cell{idx_fit,1};
    XYZ_array = current_scan(:,1:3);
    Intensity = current_scan(:,4);

    % Call the fitting function
    [C_sphere,R_sphere,E_total] = fcn_geometry_FitSphereLSQRegression(XYZ_array);


    C_sphere_array(idx_fit,:) = C_sphere;
    R_sphere_array(idx_fit,:) = R_sphere;
    E_total_array(idx_fit,:) = E_total;
    ptCloud = pointCloud(XYZ_array, 'Intensity', Intensity);
    if ~isnan(ptCloud.Location)
        view(LidarViewer,ptCloud)
        LidarViewer.Axes.Title.String = sprintf("Frame number: %d", idx_fit);
        pause(pause_time)
    end
end

%% Compute errors and visualize result
C_sphere_actual = GPS_SparkFun_Temp_ENU(fit_range,:);
% Compute the offsets between the actual sphere center and the fitting
% center
C_offset = vecnorm(C_sphere_actual-C_sphere_array,2,2);
% Compute the root mean square error (RMSE) of the center
N_fitting_center = length(fit_range);
C_rmse = sqrt(sum((C_sphere_array - C_sphere_actual).^2,2)./N_fitting_center);
figure(1221)
clf
plot(C_offset,'linewidth',2,color='blue',marker='.',markersize = 20)
ylabel('Offset of the Center [m]','FontSize',16)
figure(1222)
clf
plot(C_rmse,'linewidth',2,color='blue',marker='.',markersize = 20)
ylabel('RMSE of the Center [m]','FontSize',16)

figure(2342)
clf
fcn_Internal_showPairsOfPoints(C_sphere_actual,C_sphere_array);
legend('Temp GPS Antenna','Calibrated Target Center',fontsize=16)
xlabel('X [m]')
ylabel('Y [m]')
zlabel('Z [m]')
% axis equal


%% Internal function
function fcn_Internal_showPairsOfPoints(data_set_1,data_set_2)
    plot3(data_set_1(:,1),data_set_1(:,2),data_set_1(:,3),'linewidth',2,color='red',marker='.',markersize = 20)
    hold on
    plot3(data_set_2(:,1),data_set_2(:,2),data_set_2(:,3),'linewidth',2,color='blue',marker='.',markersize = 20)
    N_points = size(data_set_1,1);
    for idx_points = 1:N_points
        points_pair = [data_set_1(idx_points,1),data_set_1(idx_points,2),data_set_1(idx_points,3);
                        data_set_2(idx_points,1),data_set_2(idx_points,2),data_set_2(idx_points,3)];
        plot3(points_pair(:,1),points_pair(:,2),points_pair(:,3),Color='green',LineWidth=2)
    end
end