% script_test_fcn_geometry_selfCrossProduct.m
% Tests fcn_geometry_selfCrossProduct.m

% Revision history:
%      2021_04_25:
%      -- first write of the code copying functionality from 
%      fcn_FastestTraversal_checkInputsToFunctions

close all;

%% BASIC example - find the  cross-products for positive bend
fig_num = 1;
path = [0 0; 1 1; 0 2];
[cross_products] = ...
    fcn_geometry_selfCrossProduct(...
    path, fig_num) %#ok<*NOPTS,*NASGU>


%% BASIC example - find all the cross-products for negative bend
fig_num = 2;
path = [0 0; 1 1; 2 0];
[cross_products] = ...
    fcn_geometry_selfCrossProduct(...
    path, fig_num) %#ok<*NASGU>

%% BASIC example - 
% find all the cross-products for many bends
fig_num = 3;
path = [0 0; 1 1; 0 2; 2 4; 4 2; 6 2; 2 7];
[cross_products] = ...
    fcn_geometry_selfCrossProduct(...
    path, fig_num)



if 1==0 
    
    %% ERROR example - find cross product for straight
    fig_num = 4;
    path = [0 0; 1 1; 2 2];
    [cross_products,err] = ...
        fcn_geometry_selfCrossProduct(...
        path, fig_num) %#ok<*ASGLU,*NASGU>
    
    %% ERROR example - find cross product for bent back on self
    fig_num = 5;
    path = [0 0; 1 1; 0 0];
    [cross_products,err] = ...
        fcn_geometry_selfCrossProduct(...
        path, fig_num) %#ok<*NASGU>
    
    %% ERROR example - find cross product for zero-length segment
    fig_num = 6;
    path = [1 1; 1 1; 0 0];
    [cross_products,err] = ...
        fcn_geometry_selfCrossProduct(...
        path, fig_num) %#ok<*NASGU>
    
end

