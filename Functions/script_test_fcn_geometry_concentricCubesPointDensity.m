%script_test_fcn_GeomTools_concentricSquaresPointDensity.m
%This script is used to exercise the function:
%function [points] = fcn_geometry_concentricSquaresPointDensity(exterior_size,interior_size,exterior_density,interior_density,varargin)
%This function was written on 2024_6_17 by A. Goncharov, opg5041@psu.edu
%
%
%
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


 

%% BASIC example 1 - Generate points inside two cubes, no noise
%Function inputs exterior size 5, interior size 1, density of 1 point in
%exterior, 100 point density in interior.
%input variables
exterior_size=5; %units
interior_size=1; %units
exterior_density=1; %points per unit
interior_density=100; %points per unit

noise=0;
diagonal_flag=[];
fig_num=1111;

%Function
[points] = fcn_geometry_concentricSquaresPointDensity(exterior_size,interior_size,exterior_density,interior_density,noise,diagonal_flag,fig_num);


%%  BASIC example 2 - Generate points inside twocubes, noise added
%Function inputs exterior size 5, interior size 1, density of 1 point in
%exterior, 100 point density in interior, noise = 1;
%input variables
exterior_size=5; %units
interior_size=1; %units
exterior_density=1; %points per unit
interior_density=100; %points per unit

noise=1;
diagonal_flag=[];
fig_num=2222;

%Function
[points] = fcn_geometry_concentricSquaresPointDensity(exterior_size,interior_size,exterior_density,interior_density,noise,diagonal_flag,fig_num);


%% BASIC_example 3 - Create a diagonal line across, noise added
%Function inputs exterior size 5, interior size 1, density of 1 point in
%exterior, 100 point density in interior, noise = 1 and diagonal flag equal
%to 1

%input variables
exterior_size=5; %units
interior_size=1; %units
exterior_density=1; %points per unit
interior_density=100; %points per unit

noise=1;
diagonal_flag=1;
fig_num=2222;

%Function
[points] = fcn_geometry_concentricSquaresPointDensity(exterior_size,interior_size,exterior_density,interior_density,noise,diagonal_flag,fig_num);

%% BASIC_example 4 - Create a diagonal line across, noise added, fig 31
%Function inputs exterior size 5, interior size 1, density of 1 point in
%exterior, 100 point density in interior, noise = 1 and diagonal flag equal
%to 1

%input variables
exterior_size=5; %units
interior_size=1; %units
exterior_density=1; %points per unit
interior_density=100; %points per unit

noise=1;
diagonal_flag=1;
fig_num=31;

%Function
[points] = fcn_geometry_concentricSquaresPointDensity(exterior_size,interior_size,exterior_density,interior_density,noise,diagonal_flag,fig_num);

%%
