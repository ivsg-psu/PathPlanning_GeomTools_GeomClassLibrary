% script_test_fcn_geometry_findPhiConstraints
% Exercises the function: fcn_geometry_findPhiConstraints

% Revision History:
%  2020_04_06
%  -- written by Veronica Gruning (vag5076@psu.edu)
%  2020_05_02
%  -- edits by Sean Brennan to add more comments, check results,
%  add debug option
%  2021_05_27
%  -- edits for input checking, prep for geometry library

close all;


%% Simple test - should result in 90 to 0
fig_num = 1;

p_apex = [0 0];
vertex_1 = [1  0];
vertex_2 = [0 1];
[phi_start,change] = fcn_geometry_findPhiConstraints(p_apex,vertex_1,vertex_2,fig_num);
% fprintf(1,'Results of phi constrainter: \n');
% fprintf('\t %.2f \t %.2f\n', ...
%     mod(phi_start,2*pi)*180/pi, ...
%     mod(change,2*pi)*180/pi);

assert(isequal(round(phi_start,4),1.5708));
assert(isequal(round(change,4),-1.5708));

%% Simple test - different cross product direction - should result in 90 to 0
fig_num = 2;

p_apex = [0 0];
vertex_2 = [1  0];
vertex_1 = [0 1];
[phi_start,change] = fcn_geometry_findPhiConstraints(p_apex,vertex_1,vertex_2,fig_num);
% fprintf(1,'Results of phi constrainter: \n');
% fprintf('\t %.2f \t %.2f\n', ...
%     mod(phi_start,2*pi)*180/pi, ...
%     mod(change,2*pi)*180/pi);

assert(isequal(round(phi_start,4),0));
assert(isequal(round(change,4),1.5708));

%% Simple test - should result in 0 to 90
fig_num = 3;

p_apex = [0 0];
vertex_1 = [-1  0];
vertex_2 = [0 -1];
[phi_start,change] = fcn_geometry_findPhiConstraints(p_apex,vertex_1,vertex_2,fig_num);
% fprintf(1,'Results of phi constrainter: \n');
% fprintf('\t %.2f \t %.2f\n', ...
%     mod(phi_start,2*pi)*180/pi, ...
%     mod(change,2*pi)*180/pi);

assert(isequal(round(phi_start,4),4.7124));
assert(isequal(round(change,4),-1.5708));

%% Simple test - should result in 90 to 45
fig_num = 9;

p_apex = [0 0];
vertex_1 = [1  0];
vertex_2 = [-1 1];
[phi_start,change] = fcn_geometry_findPhiConstraints(p_apex,vertex_1,vertex_2,fig_num);
% fprintf(1,'Results of phi constrainter: \n');
% fprintf('\t %.2f \t %.2f\n', ...
%     mod(phi_start,2*pi)*180/pi, ...
%     mod(change,2*pi)*180/pi);

assert(isequal(round(phi_start,4),1.5708));
assert(isequal(round(change,4),-0.7854));

%% Simple test to check the vectorization - should result in 0 to 90, etc

fig_num = 11;

% -45 to 45, 180 to 270
p_apex =   [0 0;  2                      0; 0 -1; -1 0];
vertex_1 = [1  0; 2+1/(2^0.5)  -1/(2^0.5) ; -1 -1; -2 0];
vertex_2 = [0  1; 2+1/(2^0.5)   1/(2^0.5) ; 0  -2; -1 1];

[phi_start,change] = fcn_geometry_findPhiConstraints(p_apex,vertex_1,vertex_2,fig_num);
% fprintf(1,'Results of phi constrainter: \n');
% for ith_vertex = 1:length(p_apex(:,1))
%     fprintf('\t %.2f \t %.2f\n', ...
%         mod(phi_start(ith_vertex,1),2*pi)*180/pi, ...
%         mod(change(ith_vertex,1),2*pi)*180/pi);
% end

%% Simple rotation test to check for unknown items

fig_num = 13;

p_apex =   [0 0;  2                      0; 0 -1; -1 0];
vertex_1 = [1  0; 2+1/(2^0.5)  -1/(2^0.5) ; -1 -1; -2 0];
vertex_2 = [0  1; 2+1/(2^0.5)   1/(2^0.5) ; 0  -2; -1 1];

angles = 0:0.01:2*pi;
for i_angle = 1:length(angles)
    angle = angles(i_angle);
    R = [cos(angle) -sin(angle); sin(angle) cos(angle)];

    p_apex2 = p_apex*R;
    vertex_12   = vertex_1*R;
    vertex_22   = vertex_2*R;

    [phi_start,change] = fcn_geometry_findPhiConstraints(p_apex2,vertex_12,vertex_22,fig_num);
    % fprintf('\t %.2f \t %.2f\n', ...
    %     mod(phi_start,2*pi)*180/pi, ...
    %     mod(change,2*pi)*180/pi);
    drawnow;
end

%% Test of fast implementation mode 

p_apex = [0 0];
vertex_1 = [1  0];
vertex_2 = [-1 1];

% Perform the calculation in slow mode
fig_num = [];
REPS = 100; minTimeSlow = Inf; 
tic;
for i=1:REPS
    tstart = tic;
    [phi_start,change] = fcn_geometry_findPhiConstraints(p_apex,vertex_1,vertex_2, (fig_num));
    telapsed = toc(tstart);
    minTimeSlow = min(telapsed,minTimeSlow);
end
averageTimeSlow = toc/REPS;

% Perform the operation in fast mode
fig_num = -1;
minTimeFast = Inf; nsum = 10;
tic;
for i=1:REPS
    tstart = tic;
    [phi_start,change] = fcn_geometry_findPhiConstraints(p_apex,vertex_1,vertex_2, (fig_num));
    telapsed = toc(tstart);
    minTimeFast = min(telapsed,minTimeFast);
end
averageTimeFast = toc/REPS;

fprintf(1,'\n\nComparison of fast and slow modes of fcn_geometry_fitVectorToNPoints:\n');
fprintf(1,'N repetitions: %.0d\n',REPS);
fprintf(1,'Slow mode average speed per call (seconds): %.5f\n',averageTimeSlow);
fprintf(1,'Slow mode fastest speed over all calls (seconds): %.5f\n',minTimeSlow);
fprintf(1,'Fast mode average speed per call (seconds): %.5f\n',averageTimeFast);
fprintf(1,'Fast mode fastest speed over all calls (seconds): %.5f\n',minTimeFast);
fprintf(1,'Average ratio of fast mode to slow mode (unitless): %.3f\n',averageTimeSlow/averageTimeFast);
fprintf(1,'Fastest ratio of fast mode to slow mode (unitless): %.3f\n',minTimeSlow/minTimeFast);


