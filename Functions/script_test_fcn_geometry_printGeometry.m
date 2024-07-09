%% script_test_fcn_geometry_printGeometry
% Exercises the function: fcn_geometry_printGeometry

% Revision history:
% 2024_07_06 - S. Brennan
% -- wrote the code, using script_test_fcn_geometry_plotGeometry as a starter



close all;

%% BASIC test - header plotting
fid = 1;
flag_print_header = 1;
lead_string = '';

fprintf(1,'\n\n');

fcn_geometry_printGeometry('none', [], (flag_print_header), (lead_string), (fid))
fcn_geometry_printGeometry('circle', [], (flag_print_header), (lead_string), (fid))
fcn_geometry_printGeometry('arc', [], (flag_print_header), (lead_string), (fid))
fcn_geometry_printGeometry('line', [], (flag_print_header), (lead_string), (fid))
fcn_geometry_printGeometry('segment', [], (flag_print_header), (lead_string), (fid))
fcn_geometry_printGeometry('spiral', [], (flag_print_header), (lead_string), (fid))

%% BASIC test - no header plotting - default numbers
fid = 1;
flag_print_header = 0;
lead_string = 'test';

fprintf(1,'\n\n');

fcn_geometry_printGeometry('none', [], (flag_print_header), (lead_string), (fid))
fcn_geometry_printGeometry('circle', [1 2 3], (flag_print_header), (lead_string), (fid))
fcn_geometry_printGeometry('arc', [1 2 3 pi/2 3*pi/2 0 1], (flag_print_header), (lead_string), (fid))
fcn_geometry_printGeometry('line', [1 2 pi/2], (flag_print_header), (lead_string), (fid))
fcn_geometry_printGeometry('segment', [1 2 pi/2 3], (flag_print_header), (lead_string), (fid))
fcn_geometry_printGeometry('spiral', [1 2 pi/2 3 1/4 1/5], (flag_print_header), (lead_string), (fid))


%% BASIC test - no header plotting - dimensioned numbers
fid = 1;
flag_print_header = -1;
lead_string = 'test2';

fprintf(1,'\n\n');

fcn_geometry_printGeometry('none', [], (flag_print_header), (lead_string), (fid))
fcn_geometry_printGeometry('circle', [1 2 3], (flag_print_header), (lead_string), (fid))
fcn_geometry_printGeometry('arc', [1 2 3 pi/2 3*pi/2 0 1], (flag_print_header), (lead_string), (fid))
fcn_geometry_printGeometry('line', [1 2 pi/2], (flag_print_header), (lead_string), (fid))
fcn_geometry_printGeometry('segment', [1 2 pi/2 3], (flag_print_header), (lead_string), (fid))
fcn_geometry_printGeometry('spiral', [1 2 pi/2 3 1/4 1/5], (flag_print_header), (lead_string), (fid))


%% Fail conditions
if 1==0
    %% FAIL 1: points not long enough
    points = [2 3];
    [slope,intercept] = fcn_geometry_printGeometry(points,fid);
    fprintf(1,'\n\nSlope is: %.2f, Intercept is: %.2f\n',slope,intercept);
end