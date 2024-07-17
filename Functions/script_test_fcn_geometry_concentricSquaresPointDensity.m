%script_test_fcn_geometry_concentricSquaresPointDensity.m
%This script is used to exercise the function:
%fcn_geometry_concentricSquaresPointDensity
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
Int_Square_Size=3;
fig_num = 1111;

[points] = fcn_geometry_concentricSquaresPointDensity(N_points,Ext_Square_Size,Int_Square_Size,[],[],fig_num);

assert(length(points)==N_points);


%%  BASIC example 2 - Generate points inside two squares, noise added

N_points = 100;
Ext_Square_Size=10;
Int_Square_Size=3;
fig_num = 2222;
noise= 0.5;

[points] = fcn_geometry_concentricSquaresPointDensity(N_points,Ext_Square_Size,Int_Square_Size,noise,[],fig_num);
assert(length(points)==N_points);


%% BASIC_example 3 - Create a diagonal line across, no noise

N_points = 100;
Ext_Square_Size=10;
Int_Square_Size=3;
fig_num = 3333;
diag_flag=1;

[points] = fcn_geometry_concentricSquaresPointDensity(N_points,Ext_Square_Size,Int_Square_Size,[],diag_flag,fig_num);
assert(length(points)==N_points);

%% BASIC_example 4 - Create a diagonal line across, noise added

N_points = 100;
Ext_Square_Size=10;
Int_Square_Size=3;
fig_num = 4444;
diag_flag=1;
noise= 0.3;

[points] = fcn_geometry_concentricSquaresPointDensity(N_points,Ext_Square_Size,Int_Square_Size,noise,diag_flag,fig_num);
assert(length(points)==N_points);


%% BASIC_example 5 - Create a diagonal line across, noise added, fig 31


N_points = 100;
Ext_Square_Size=10;
Int_Square_Size=3;
fig_num = 31;
diag_flag=1;
noise= 0.3;

[points] = fcn_geometry_concentricSquaresPointDensity(N_points,Ext_Square_Size,Int_Square_Size,noise,diag_flag,fig_num);
assert(length(points)==N_points);


%%
