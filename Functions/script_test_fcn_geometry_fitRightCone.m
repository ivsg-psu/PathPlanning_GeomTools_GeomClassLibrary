
% script_test_fcn_geometry_fitRightCone.m
% This is a script to exercise the function: fcn_geometry_fitRightCone.m
% This function was written on 2023_01_22 by Xinyu Cao, xfc5113@psu.edu

% Revision history:
% 2024_01_24
% -- first write of the code

%% Set up the workspace
close all

%% Load Data: Use the rawdata_lidar.mat in the Data folder to start, 
% the loaded data is a struct containing two fields
% rawdata.GPS_SparkFun_Temp_ENU and rawdata.Lidar_pointcloud_cell
load('PointCloud_Separated_Data_newTarget.mat')

% Test 1: a basic test, change the index of scan and ring label if needed
idx_scan = 1;
test_ring = ptCloud_pts_layers_separated_cell{idx_scan}.Ring0;
inputPoints = test_ring(:,1:3);
ring_id = test_ring(1,5);
fig_num = 1;
[cone_parameters,fittedPoints,fitting_result] = fcn_geometry_fitRightCone(inputPoints,ring_id,fig_num);
