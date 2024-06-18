%script_test_fcn_GeomTools_concentricSquaresPointDensity.m
%This script is used to exercise the function:
%fcn_GeomTools_concentricSquaresPointsDenisty
%This function was written on 2024_6_17 by A. Goncharov, opg5041@psu.edu
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
%% BASIC example 1 - Generate points inside two squares, no noise
%Function inputs 100 points, and an exterior square size 10.

N_points = 100;
Ext_Square_Size=10;
fig_num = 1111;

[X,Y,Z] = fcn_geometry_concentricSquaresPointDensity(N_points,Ext_Square_Size,[],[],fig_num);

assert(length(X)==100);


%%  BASIC example 2 - Generate points inside two squares, noise added

fcn_geometry_concentricSquaresPointDensity(100,10,0.5)


%% BASIC_example 3 - Create a diagonal line across, no noise

fcn_geometry_concentricSquaresPointDensity(100,10,0,1)


%% BASIC_example 4 - Create a diagonal line across, noise added

fcn_geometry_concentricSquaresPointDensity(100,10,0.5,1)


%% BASIC_example 5 - Create a diagonal line across, noise added, fig 31


fcn_geometry_concentricSquaresPointDensity(100,10,[],[],31)


%%
