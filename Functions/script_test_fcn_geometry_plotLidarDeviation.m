%script_test_fcn_geometry_plotLidarDeviation.m
%This script is used to exercise the function:
%fcn_geometry_plotLidarDeviation
%This function was written on 2024_7_17 by A. Goncharov, opg5041@psu.edu
%
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

%% Basic Example 1 - Load Data and Generate Three Distinct Figures Store Variables

cell_start=10;
cell_end=20;
color_map='jet';
data_name='TestTrack_Entire_Loop';
fig_1=1111;
fig_2=2222;
fig_3=3333;

[z,mean_z,deviation_from_mean,angles,mean_angle,angle_difference] ...
= fcn_geometry_plotLidarDeviation(cell_start,cell_end,color_map,...
data_name,fig_1,fig_2,fig_3);

%% Basic Example 2 - Load Data and Generate Three Distinct Figures Don't Store

cell_start=10;
cell_end=20;
color_map='jet';
data_name='TestTrack_Entire_Loop';
fig_1=[];
fig_2=[];
fig_3=3234;

fcn_geometry_plotLidarDeviation(cell_start,cell_end,color_map,...
data_name,fig_1,fig_2,fig_3);

%% Basic Example 3 - Load Data, No Figures

cell_start=10;
cell_end=20;
color_map='jet';
data_name='TestTrack_Entire_Loop';
fig_1=[];
fig_2=[];
fig_3=[];

fcn_geometry_plotLidarDeviation(cell_start,cell_end,color_map,...
data_name,fig_1,fig_2,fig_3);

%% Basic Example 4 - Load Data, Two Figures, Different Color Map, Different cell_start

cell_start=1;
cell_end=30;
color_map='prism';
data_name='TestTrack_Entire_Loop';
fig_1=1111;
fig_2=2222;
fig_3=[];

fcn_geometry_plotLidarDeviation(cell_start,cell_end,color_map,...
data_name,fig_1,fig_2,fig_3);