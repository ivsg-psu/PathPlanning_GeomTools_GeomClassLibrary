%Script for fcn_geometry_concatenatePoints
%Written by Aleksandr Goncharov
%Revision History:
% -- 7/29/2024 -- Aleksandr Goncharov
% -- Updated examples to include reference LLA coordinates

% [concatenate_LiDAR_XYZ_points,concatenate_VehiclePose_XYZ_points...
%     ,concatenate_unit_ortho_vehicle_vectors_XYZ] = ...
%     fcn_geometry_concatenatePoints(LiDAR_Scan_ENU_Entire_Loop,VehiclePose,...
%     rings_to_analyze,scanLineNumber_start,scanLineNumber_end,unit_ortho_vehicle_vectors_XY,varargin)

%Load the unit_ortho_vehicle_vector
Color_Track=[1 1 0];
Color_Start=[0 1 0];
Color_End=[1 0 0];
Marker_Size_Track= 10;
Marker_Size_Ends= 10;

[vehicle_change_in_pose_XY,unit_vehicle_change_in_pose_XY,unit_ortho_vehicle_vectors_XY] = ...
    fcn_geometry_vehiclePosition(VehiclePose,Color_Track,Color_Start,Color_End,Marker_Size_Track,Marker_Size_Ends);


%% Basic Example
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   ____            _        ______                           _
%  |  _ \          (_)      |  ____|                         | |
%  | |_) | __ _ ___ _  ___  | |__  __  ____ _ _ __ ___  _ __ | | ___
%  |  _ < / _` / __| |/ __| |  __| \ \/ / _` | '_ ` _ \| '_ \| |/ _ \
%  | |_) | (_| \__ \ | (__  | |____ >  < (_| | | | | | | |_) | |  __/
%  |____/ \__,_|___/_|\___| |______/_/\_\__,_|_| |_| |_| .__/|_|\___|
%                                                      | |
%                                                      |_|
% See: https://patorjk.com/software/taag/#p=display&f=Big&t=Basic%20Example
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%§

%% Example 1: Rings 1-3, lines 1400-1410

rings_to_analyze = [1 2 3];
scanLineNumber_start=1400;
scanLineNumber_end = 1410;
figure_num=1111;
ENU_3D_fig_num=3;
simple_view_toggle_1=[];
scaling_1=3;
ENU_XZ_fig_num=4;
ENU_YZ_fig_num=5;
scaling_2=3;
reference_latitude = 40.86368573;
reference_longitude = -77.83592832;
reference_altitude = 344.189;
LLA_fig_num=2;

[concatenate_LiDAR_XYZ_points,concatenate_VehiclePose_XYZ_points...
    ,concatenate_unit_ortho_vehicle_vectors_XYZ,intensity_min, intensity_max] = fcn_geometry_concatenatePoints(LiDAR_Scan_ENU_Entire_Loop,VehiclePose,...
    rings_to_analyze,scanLineNumber_start,scanLineNumber_end,unit_ortho_vehicle_vectors_XY,figure_num,ENU_3D_fig_num,simple_view_toggle_1,scaling_1,...
    ENU_XZ_fig_num,ENU_YZ_fig_num,scaling_2,reference_latitude,reference_longitude,reference_altitude,LLA_fig_num);
assert(length(concatenate_unit_ortho_vehicle_vectors_XYZ)==length(concatenate_VehiclePose_XYZ_points))

%% Example 2: Ring 4, lines 1400-1410

rings_to_analyze = 4;
scanLineNumber_start=1400;
scanLineNumber_end = 1410;
figure_num=2222;
ENU_3D_fig_num=3;
simple_view_toggle_1=[];
scaling_1=3;
ENU_XZ_fig_num=4;
ENU_YZ_fig_num=5;
scaling_2=3;
reference_latitude = 40.86368573;
reference_longitude = -77.83592832;
reference_altitude = 344.189;
LLA_fig_num=2;

[concatenate_LiDAR_XYZ_points,concatenate_VehiclePose_XYZ_points...
    ,concatenate_unit_ortho_vehicle_vectors_XYZ,intensity_min, intensity_max] = fcn_geometry_concatenatePoints(LiDAR_Scan_ENU_Entire_Loop,VehiclePose,...
    rings_to_analyze,scanLineNumber_start,scanLineNumber_end,unit_ortho_vehicle_vectors_XY,figure_num,ENU_3D_fig_num,simple_view_toggle_1,scaling_1,...
    ENU_XZ_fig_num,ENU_YZ_fig_num,scaling_2,reference_latitude,reference_longitude,reference_altitude,LLA_fig_num);
assert(length(concatenate_unit_ortho_vehicle_vectors_XYZ)==length(concatenate_VehiclePose_XYZ_points))

%% Example 3: Ring 4, lines 1200-1300
rings_to_analyze = 4;
scanLineNumber_start=1200;
scanLineNumber_end = 1300;
figure_num=3333;
ENU_3D_fig_num=3;
simple_view_toggle_1=[];
scaling_1=3;
ENU_XZ_fig_num=4;
ENU_YZ_fig_num=5;
scaling_2=3;
reference_latitude = 40.86368573;
reference_longitude = -77.83592832;
reference_altitude = 344.189;
LLA_fig_num=2;

[concatenate_LiDAR_XYZ_points,concatenate_VehiclePose_XYZ_points...
    ,concatenate_unit_ortho_vehicle_vectors_XYZ,intensity_min, intensity_max] = fcn_geometry_concatenatePoints(LiDAR_Scan_ENU_Entire_Loop,VehiclePose,...
    rings_to_analyze,scanLineNumber_start,scanLineNumber_end,unit_ortho_vehicle_vectors_XY,figure_num,ENU_3D_fig_num,simple_view_toggle_1,scaling_1,...
    ENU_XZ_fig_num,ENU_YZ_fig_num,scaling_2,reference_latitude,reference_longitude,reference_altitude,LLA_fig_num);
assert(length(concatenate_unit_ortho_vehicle_vectors_XYZ)==length(concatenate_VehiclePose_XYZ_points))


