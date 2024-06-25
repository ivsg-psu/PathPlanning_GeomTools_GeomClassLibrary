%script_test_fcn_GeomTools_concentricSquaresPointDensity.m
%This script is used to exercise the function:
%fcn_geometry_concentricSquaresPointsDenisty
%This function was written on 2024_6_17 by A. Goncharov, opg5041@psu.edu
%
%Revision History:
% 2024_6_24 - A. Goncharov - Revised to have point concentration input instead of N_points
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

ext_length=10;
int_length=5;
ext_point_concentration=3;
int_point_concentration=20;

fig_num = 1111;

[points] = fcn_geometry_concentricSquaresPointDensity(ext_length,int_length,ext_point_concentration,int_point_concentration,[],[],fig_num);

assert(length(points)==pointChecker(ext_length,int_length,ext_point_concentration,int_point_concentration));


%%  BASIC example 2 - Generate points inside two squares, noise added

ext_length=10;
int_length=5;
ext_point_concentration=3;
int_point_concentration=20;
noise=1;

fig_num = 2222;

[points] = fcn_geometry_concentricSquaresPointDensity(ext_length,int_length,ext_point_concentration,int_point_concentration,noise,[],fig_num);

assert(length(points)==pointChecker(ext_length,int_length,ext_point_concentration,int_point_concentration));



%% BASIC_example 3 - Create a diagonal line across, no noise

ext_length=10;
int_length=5;
ext_point_concentration=3;
int_point_concentration=20;
diagonal_flag=1;

fig_num = 3333;

[points] = fcn_geometry_concentricSquaresPointDensity(ext_length,int_length,ext_point_concentration,int_point_concentration,[],diagonal_flag,fig_num);

assert(length(points)==pointChecker(ext_length,int_length,ext_point_concentration,int_point_concentration));

%% BASIC_example 4 - Create a diagonal line across, noise added

ext_length=10;
int_length=5;
ext_point_concentration=3;
int_point_concentration=20;
noise=1;
diagonal_flag=1;
fig_num = 4444;

[points] = fcn_geometry_concentricSquaresPointDensity(ext_length,int_length,ext_point_concentration,int_point_concentration,noise,diagonal_flag,fig_num);

assert(length(points)==pointChecker(ext_length,int_length,ext_point_concentration,int_point_concentration));

%% BASIC_example 5 - Create a diagonal line across, noise added, fig 31


ext_length=10;
int_length=5;
ext_point_concentration=3;
int_point_concentration=20;
noise=1;
diagonal_flag=1;
fig_num = 31;

[points] = fcn_geometry_concentricSquaresPointDensity(ext_length,int_length,ext_point_concentration,int_point_concentration,noise,diagonal_flag,fig_num);

assert(length(points)==pointChecker(ext_length,int_length,ext_point_concentration,int_point_concentration));


%%





function [numPoints] = pointChecker(ext_length,int_length,ext_point_concentration,int_point_concentration)

ext_points= round(ext_point_concentration*ext_length^2);
int_points= round(int_point_concentration*int_length^2);

n_ext = round(sqrt(ext_points));
n_int = round(sqrt(int_points));

ext_side=ext_length/2;
int_side=int_length/2;

%grid generation
[x_e,y_e]=meshgrid(linspace(-ext_side,ext_side,n_ext),linspace(-ext_side,ext_side,n_ext));
x_ext=x_e(:);
y_ext=y_e(:);

[x_i,y_i]=meshgrid(linspace(-int_side,int_side,n_int),linspace(-int_side,int_side,n_int),0);
x_int=x_i(:);
y_int=y_i(:);


%point remover
interior_mask= abs(x_ext) <= int_side & abs(y_ext) <= int_side;

x_ext=x_ext(~interior_mask);
y_ext=y_ext(~interior_mask);
z_ext=zeros(length(x_ext),1);
z_int=zeros(length(x_int),1);

X=[x_ext;x_int];
Y=[y_ext;y_int];
Z=[z_ext;z_int];

points=[X Y Z];

numPoints=length(points);

end