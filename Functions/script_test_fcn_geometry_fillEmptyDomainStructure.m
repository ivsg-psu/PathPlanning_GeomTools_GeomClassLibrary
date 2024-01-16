% script_test_fcn_geometry_fillEmptyDomainStructure
% Exercises the function: fcn_geometry_fillEmptyDomainStructure
% Revision history:
% 2024_01_15
% -- wrote the code


%% Test 1: a basic test of line segment fitting, specifying index-type base_point_index
emptyDomain = fcn_geometry_fillEmptyDomainStructure;


%% Fail conditions
if 1==0
    %% FAIL 1: points not long enough
    points = [2 3];
    [slope,intercept] = fcn_geometry_fillEmptyDomainStructure(points,fig_num);
    fprintf(1,'\n\nSlope is: %.2f, Intercept is: %.2f\n',slope,intercept);
end