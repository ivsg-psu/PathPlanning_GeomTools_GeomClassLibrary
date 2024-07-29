%% script_test_fcn_geometry_vehicleOrientation
% Script for function: fcn_geometry_vehicleOrientation
% Written on: 7_26_2024
% By: Aleksandr Goncharov - opg5041@psu.edu
% 
% Revision History:

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
%% Example 1: All features utilized
Color_Track=[1 1 0];
Color_Start=[0 1 0];
Color_End=[1 0 0];
Marker_Size_Track= 10;
Marker_Size_Ends= 10;
Fig_Num=1234;

[test1_vehicle,test1_unit,test1_ortho] = ...
    fcn_geometry_vehiclePosition(VehiclePose,Color_Track,Color_Start,Color_End,Marker_Size_Track,Marker_Size_Ends,Fig_Num);
assert(length(test1_vehicle)==length(VehiclePose))


%% Example 2: Only the default features and figure number

Color_Track=[];
Color_Start=[];
Color_End=[];
Marker_Size_Track=[];
Marker_Size_Ends=[];
Fig_Num=2222;

[test2_vehicle,test2_unit,test2_ortho]=fcn_geometry_vehiclePosition(VehiclePose,Color_Track,Color_Start,Color_End,Marker_Size_Track,Marker_Size_Ends,Fig_Num);
assert(length(vehicle_change_in_pose_XY)==length(VehiclePose))

assert(length(test2_vehicle)==length(VehiclePose))

%% Example 3: Only the default features and color changed for ends and figure number

Color_Track=[];
Color_Start=[];
Color_End=[1 0 1];
Marker_Size_Track=[];
Marker_Size_Ends=[];
Fig_Num=3333;

[test3_vehicle,test3_unit,test3_ortho]=fcn_geometry_vehiclePosition(VehiclePose,Color_Track,Color_Start,Color_End,Marker_Size_Track,Marker_Size_Ends,Fig_Num);
assert(length(vehicle_change_in_pose_XY)==length(VehiclePose))

assert(length(test3_vehicle)==length(VehiclePose))

%% Example 4: Big marker sizes

Color_Track=[];
Color_Start=[];
Color_End=[];
Marker_Size_Track=10;
Marker_Size_Ends=50;
Fig_Num=4444;

[test4_vehicle,test4_unit,test4_ortho]=fcn_geometry_vehiclePosition(VehiclePose,Color_Track,Color_Start,Color_End,Marker_Size_Track,Marker_Size_Ends,Fig_Num);
assert(length(vehicle_change_in_pose_XY)==length(VehiclePose))

assert(length(test4_vehicle)==length(VehiclePose))


