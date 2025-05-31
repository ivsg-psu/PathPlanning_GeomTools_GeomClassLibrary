% script_test_fcn_geometry_findPlaneNormal
% Exercises the function: fcn_geometry_findPlaneNormal
% Revision history:
% 2025_05_30
% -- wrote the code

close all;


%% Demonstration Examples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  _____                                 _             _   _               ______                           _
% |  __ \                               | |           | | (_)             |  ____|                         | |
% | |  | | ___ _ __ ___   ___  _ __  ___| |_ _ __ __ _| |_ _  ___  _ __   | |__  __  ____ _ _ __ ___  _ __ | | ___  ___
% | |  | |/ _ \ '_ ` _ \ / _ \| '_ \/ __| __| '__/ _` | __| |/ _ \| '_ \  |  __| \ \/ / _` | '_ ` _ \| '_ \| |/ _ \/ __|
% | |__| |  __/ | | | | | (_) | | | \__ \ |_| | | (_| | |_| | (_) | | | | | |____ >  < (_| | | | | | | |_) | |  __/\__ \
% |_____/ \___|_| |_| |_|\___/|_| |_|___/\__|_|  \__,_|\__|_|\___/|_| |_| |______/_/\_\__,_|_| |_| |_| .__/|_|\___||___/
%                                                                                                    | |
%                                                                                                    |_|
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Demonstration%20Examples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Demonstration case 1: basic plane
fig_num = 0001;
figure(fig_num);
clf;

points = [0 0 0; 0 1 0; 1 0 0];

[unit_normal_vector, base_point, flags_in_agreement] = fcn_geometry_findPlaneNormal(points,(fig_num));

% Check variable types
assert(isnumeric(unit_normal_vector));
assert(isnumeric(base_point));
assert(isnumeric(flags_in_agreement));


% Check variable lengths
Npoints = length(points(:,1));
assert(isequal(size(parameters),[1 4]));
assert(isequal(size(standard_deviation_in_z),[1 1]));
assert(isequal(size(z_fit),[Npoints 1]));
assert(isequal(size(unit_vector),[1 3]));
assert(isequal(size(standard_deviation_in_plane_orthogonals),[1 1]));

% Check variable values
assert(isequal(round(parameters,4),round(true_parameters,4)));
assert(isequal(round(standard_deviation_in_z,4),0));
assert(isequal(round(z_fit,4),round(true_z,4)));
assert(isequal(round(unit_vector,4),round(true_parameters(1,1:3),4)));
assert(isequal(round(standard_deviation_in_plane_orthogonals,4),0));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Demonstration case 2: basic plane with noise
fig_num = 0002;
figure(fig_num);
clf;


Npoints = 100;

unNormal_true_parameters = [.1 .2 1 3];
vectorLength = sum(unNormal_true_parameters(1,1:3).^2,2).^0.5;
true_parameters = unNormal_true_parameters./vectorLength;

Constant = true_parameters(4);
% X = [0; 1; 0];
% Y = [0; 1; 1];
% Z = (Constant*ones(3,1)  - X*true_parameters(1,1) - Y*true_parameters(1,2))./true_parameters(1,3);
% points = [X Y Z];


points = randn(Npoints,3);
X = points(:,1);
Y = points(:,2);
true_z = (Constant*ones(Npoints,1)  - X*true_parameters(1,1) - Y*true_parameters(1,2))./true_parameters(1,3); % Solve for z vertices data

true_sigma = 0.001;
Z = true_z + true_sigma*randn(Npoints,1);

[fitted_parameters, standard_deviation_in_z] = fcn_geometry_findPlaneNormal([X Y Z],fig_num);

assert(isequal(round(true_parameters,1),round(fitted_parameters,1)));
assert(isequal(round(true_sigma,1),round(standard_deviation_in_z,1)));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%%
% From: https://www.mathworks.com/matlabcentral/fileexchange/43305-plane-fit
fig_num = 0003;
figure(fig_num);
clf;

% Generate points that lie approximately in the Z=0 plane
N = 10;
[X,Y] = meshgrid(linspace(0,1,N));
XYZ_1 = [X(:) Y(:) 0.05*randn(N^2,1)];
plot3(XYZ_1(:,1),XYZ_1(:,2),XYZ_1(:,3),'r.');
hold on

%compute the normal to the plane and a point that belongs to the plane
% [n_1_old,~,p_1_old] = affine_fit(XYZ_1);
[~, ~, ~, n_1, p_1, ~, ~, V_1] = fcn_geometry_findPlaneNormal(XYZ_1);
n_1 = n_1'; % Original code uses transpose form


%generate points that lie approximately in the Z=X plane
%the normal vector is
n_2_exact = [-sqrt(2)/2 0 sqrt(2)/2];
N = 12;
[X,Y] = meshgrid(linspace(0,1,N));
XYZ_2 = [X(:) Y(:) X(:)] + bsxfun(@times,0.05*randn(N^2,1),n_2_exact);
plot3(XYZ_2(:,1),XYZ_2(:,2),XYZ_2(:,3),'b.');


%compute the normal to the plane and a point that belongs to the plane
% [n_2_old,V_2_old,p_2_old] = affine_fit(XYZ_2);
[~, ~, ~, n_2, p_2, ~, ~, V_2] = fcn_geometry_findPlaneNormal(XYZ_2);
n_2 = n_2'; % Original code uses transpose form

%plot the two points p_1 and p_2
plot3(p_1(1),p_1(2),p_1(3),'ro','markersize',15,'markerfacecolor','red');
plot3(p_2(1),p_2(2),p_2(3),'bo','markersize',15,'markerfacecolor','blue');

%plot the normal vector
quiver3(p_1(1),p_1(2),p_1(3),n_1(1)/3,n_1(2)/3,n_1(3)/3,'r','linewidth',2)
quiver3(p_2(1),p_2(2),p_2(3),n_2(1)/3,n_2(2)/3,n_2(3)/3,'b','linewidth',2)

%plot the two adjusted planes
[X,Y] = meshgrid(linspace(0,1,3));

%first plane
surf(X,Y, - (n_1(1)/n_1(3)*X+n_1(2)/n_1(3)*Y-dot(n_1,p_1)/n_1(3)),'facecolor','red','facealpha',0.5);

%second plane
%NB: if the plane is vertical the above method cannot be used, one should
%use the secont output of affine_fit which contains a base of the plane.
%this is illustrated below
%S1 and S2 are the coordinates of the plane points in the basis made of the
%columns ov V_2
[S1,S2] = meshgrid([-1 0 1]);
%generate the pont coordinates
X = p_2(1)+[S1(:) S2(:)]*V_2(1,:)';
Y = p_2(2)+[S1(:) S2(:)]*V_2(2,:)';
Z = p_2(3)+[S1(:) S2(:)]*V_2(3,:)';
%plot the plane
surf(reshape(X,3,3),reshape(Y,3,3),reshape(Z,3,3),'facecolor','blue','facealpha',0.5);

xlabel('x');
ylabel('y');
zlabel('z');
axis equal
%compute the angle between the planes in [0 90] degrees
angle = acosd(dot(n_1,n_2));
if angle>90
    angle = 180-angle;
end
disp(angle);




%% Basic testing examples in 3D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  ____            _        _______        _   _                ______                           _                       ____  _____
% |  _ \          (_)      |__   __|      | | (_)              |  ____|                         | |                     |___ \|  __ \
% | |_) | __ _ ___ _  ___     | | ___  ___| |_ _ _ __   __ _   | |__  __  ____ _ _ __ ___  _ __ | | ___  ___    ______    __) | |  | |
% |  _ < / _` / __| |/ __|    | |/ _ \/ __| __| | '_ \ / _` |  |  __| \ \/ / _` | '_ ` _ \| '_ \| |/ _ \/ __|  |______|  |__ <| |  | |
% | |_) | (_| \__ \ | (__     | |  __/\__ \ |_| | | | | (_| |  | |____ >  < (_| | | | | | | |_) | |  __/\__ \            ___) | |__| |
% |____/ \__,_|___/_|\___|    |_|\___||___/\__|_|_| |_|\__, |  |______/_/\_\__,_|_| |_| |_| .__/|_|\___||___/           |____/|_____/
%                                                       __/ |                             | |
%                                                      |___/                              |_|
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Basic%20Testing%20%20Examples%20%20-%203D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Basic test case 1: basic plane with 3 points, XY plane
fig_num = 1001;
figure(fig_num);
clf;

points = [0 0 0; 1 0 0; 0 1 0];

[parameters, standard_deviation_in_z, z_fit, unit_vector, base_point, standard_deviation_in_plane_orthogonals] = fcn_geometry_findPlaneNormal(points,fig_num);

% Check variable types
assert(isnumeric(parameters));
assert(isnumeric(standard_deviation_in_z));
assert(isnumeric(z_fit));
assert(isnumeric(unit_vector));
assert(isnumeric(base_point));
assert(isnumeric(standard_deviation_in_plane_orthogonals));

% Check variable lengths
Npoints = length(points(:,1));
assert(isequal(size(parameters),[1 4]));
assert(isequal(size(standard_deviation_in_z),[1 1]));
assert(isequal(size(z_fit),[Npoints 1]));
assert(isequal(size(unit_vector),[1 3]));
assert(isequal(size(standard_deviation_in_plane_orthogonals),[1 1]));

% Check variable values
assert(isequal(round(parameters,4),[0 0 1 0]));
assert(isequal(round(standard_deviation_in_z,4),0));
assert(isequal(round(z_fit,4),zeros(Npoints,1)));
assert(isequal(round(unit_vector,4),[0 0 1]));
assert(isequal(round(standard_deviation_in_plane_orthogonals,4),0));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Basic test case 2: basic plane with 4 points, XY plane
fig_num = 1002;
figure(fig_num);
clf;

points = [0 0 0; 1 0 0; 0 1 0; 2 2 0];

[parameters, standard_deviation_in_z, z_fit, unit_vector, base_point, standard_deviation_in_plane_orthogonals] = fcn_geometry_findPlaneNormal(points,fig_num);

% Check variable types
assert(isnumeric(parameters));
assert(isnumeric(standard_deviation_in_z));
assert(isnumeric(z_fit));
assert(isnumeric(unit_vector));
assert(isnumeric(base_point));
assert(isnumeric(standard_deviation_in_plane_orthogonals));

% Check variable lengths
Npoints = length(points(:,1));
assert(isequal(size(parameters),[1 4]));
assert(isequal(size(standard_deviation_in_z),[1 1]));
assert(isequal(size(z_fit),[Npoints 1]));
assert(isequal(size(unit_vector),[1 3]));
assert(isequal(size(standard_deviation_in_plane_orthogonals),[1 1]));

% Check variable values
assert(isequal(round(parameters,4),[0 0 1 0]));
assert(isequal(round(standard_deviation_in_z,4),0));
assert(isequal(round(z_fit,4),zeros(Npoints,1)));
assert(isequal(round(unit_vector,4),[0 0 1]));
assert(isequal(round(standard_deviation_in_plane_orthogonals,4),0));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Basic test case 3: basic plane with 3 points, XZ plane
fig_num = 1003;
figure(fig_num);
clf;

points = [0 0 0; 1 0 0; 0 0 1];

[parameters, standard_deviation_in_z, z_fit, unit_vector, base_point, standard_deviation_in_plane_orthogonals] = fcn_geometry_findPlaneNormal(points,fig_num);

% Check variable types
assert(isnumeric(parameters));
assert(isnumeric(standard_deviation_in_z));
assert(isnumeric(z_fit));
assert(isnumeric(unit_vector));
assert(isnumeric(base_point));
assert(isnumeric(standard_deviation_in_plane_orthogonals));

% Check variable lengths
Npoints = length(points(:,1));
assert(isequal(size(parameters),[1 4]));
assert(isequal(size(standard_deviation_in_z),[1 1]));
assert(isequal(size(z_fit),[Npoints 1]));
assert(isequal(size(unit_vector),[1 3]));
assert(isequal(size(standard_deviation_in_plane_orthogonals),[1 1]));

% Check variable values
assert(isequal(round(parameters,4),[0 1 0 0]));
assert(isnan(round(standard_deviation_in_z,4)));
assert(all(isnan((round(z_fit,4)))));
assert(isequal(round(unit_vector,4),[0 1 0]));
assert(isequal(round(standard_deviation_in_plane_orthogonals,4),0));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));


%% Basic test case 4: basic plane with 3 points, YZ plane
fig_num = 1004;
figure(fig_num);
clf;

points = [0 0 0; 0 1 0; 0 0 1];

[parameters, standard_deviation_in_z, z_fit, unit_vector, base_point, standard_deviation_in_plane_orthogonals] = fcn_geometry_findPlaneNormal(points,fig_num);

% Check variable types
assert(isnumeric(parameters));
assert(isnumeric(standard_deviation_in_z));
assert(isnumeric(z_fit));
assert(isnumeric(unit_vector));
assert(isnumeric(base_point));
assert(isnumeric(standard_deviation_in_plane_orthogonals));

% Check variable lengths
Npoints = length(points(:,1));
assert(isequal(size(parameters),[1 4]));
assert(isequal(size(standard_deviation_in_z),[1 1]));
assert(isequal(size(z_fit),[Npoints 1]));
assert(isequal(size(unit_vector),[1 3]));
assert(isequal(size(standard_deviation_in_plane_orthogonals),[1 1]));

% Check variable values
assert(isequal(round(parameters,4),[1 0 0 0]));
assert(isnan(round(standard_deviation_in_z,4)));
assert(all(isnan((round(z_fit,4)))));
assert(isequal(round(unit_vector,4),[1 0 0]));
assert(isequal(round(standard_deviation_in_plane_orthogonals,4),0));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));


%% Basic test case 5: basic plane with 3 points, X=Y plane with offset
fig_num = 1005;
figure(fig_num);
clf;

% z = (-Ax + -By + D)/Constant

unNormal_true_parameters = [.1 .2 1 4];
vectorLength = sum(unNormal_true_parameters(1,1:3).^2,2).^0.5;
true_parameters = unNormal_true_parameters./vectorLength;

Constant = true_parameters(4);
X = [0; 1; 0];
Y = [0; 1; 1];
Z = (Constant*ones(3,1)  - X*true_parameters(1,1) - Y*true_parameters(1,2))./true_parameters(1,3);

points = [X Y Z];

% Make sure equation adds up to constant: sum(points.*true_parameters(1,1:3),2)

[parameters, standard_deviation_in_z, z_fit, unit_vector, base_point, standard_deviation_in_plane_orthogonals] = fcn_geometry_findPlaneNormal(points,fig_num);

% Check variable types
assert(isnumeric(parameters));
assert(isnumeric(standard_deviation_in_z));
assert(isnumeric(z_fit));
assert(isnumeric(unit_vector));
assert(isnumeric(base_point));
assert(isnumeric(standard_deviation_in_plane_orthogonals));

% Check variable lengths
Npoints = length(points(:,1));
assert(isequal(size(parameters),[1 4]));
assert(isequal(size(standard_deviation_in_z),[1 1]));
assert(isequal(size(z_fit),[Npoints 1]));
assert(isequal(size(unit_vector),[1 3]));
assert(isequal(size(standard_deviation_in_plane_orthogonals),[1 1]));

% Check variable values
assert(isequal(round(parameters,4),round(true_parameters,4)));
assert(isequal(round(standard_deviation_in_z,4),0));
assert(isequal(round(z_fit,4),round(Z,4)));
assert(isequal(round(unit_vector,4),round(true_parameters(1,1:3),4)));
assert(isequal(round(standard_deviation_in_plane_orthogonals,4),0));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Fast Mode Tests
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  ______        _     __  __           _        _______        _
% |  ____|      | |   |  \/  |         | |      |__   __|      | |
% | |__ __ _ ___| |_  | \  / | ___   __| | ___     | | ___  ___| |_ ___
% |  __/ _` / __| __| | |\/| |/ _ \ / _` |/ _ \    | |/ _ \/ __| __/ __|
% | | | (_| \__ \ |_  | |  | | (_) | (_| |  __/    | |  __/\__ \ |_\__ \
% |_|  \__,_|___/\__| |_|  |_|\___/ \__,_|\___|    |_|\___||___/\__|___/
%
%
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Fast%20Mode%20Tests
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Basic example - NO FIGURE
fig_num = 9901;
figure(fig_num);
close(fig_num);

points = [0 0 0; 1 0 0; 0 1 0];

[parameters, standard_deviation_in_z, z_fit, unit_vector, base_point, standard_deviation_in_plane_orthogonals] = fcn_geometry_findPlaneNormal(points,[]);

% Check variable types
assert(isnumeric(parameters));
assert(isnumeric(standard_deviation_in_z));
assert(isnumeric(z_fit));
assert(isnumeric(unit_vector));
assert(isnumeric(base_point));
assert(isnumeric(standard_deviation_in_plane_orthogonals));

% Check variable lengths
Npoints = length(points(:,1));
assert(isequal(size(parameters),[1 4]));
assert(isequal(size(standard_deviation_in_z),[1 1]));
assert(isequal(size(z_fit),[Npoints 1]));
assert(isequal(size(unit_vector),[1 3]));
assert(isequal(size(standard_deviation_in_plane_orthogonals),[1 1]));

% Check variable values
assert(isequal(round(parameters,4),[0 0 1 0]));
assert(isequal(round(standard_deviation_in_z,4),0));
assert(isequal(round(z_fit,4),zeros(Npoints,1)));
assert(isequal(round(unit_vector,4),[0 0 1]));
assert(isequal(round(standard_deviation_in_plane_orthogonals,4),0));

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));

%% Basic example of vertex calculation - non-normal wall shrinking, NO FIGURE, FAST MODE
fig_num = 9902;
figure(fig_num);
close(fig_num);

points = [0 0 0; 1 0 0; 0 1 0];

[parameters, standard_deviation_in_z, z_fit, unit_vector, base_point, standard_deviation_in_plane_orthogonals] = fcn_geometry_findPlaneNormal(points,-1);

% Check variable types
assert(isnumeric(parameters));
assert(isnumeric(standard_deviation_in_z));
assert(isnumeric(z_fit));
assert(isnumeric(unit_vector));
assert(isnumeric(base_point));
assert(isnumeric(standard_deviation_in_plane_orthogonals));

% Check variable lengths
Npoints = length(points(:,1));
assert(isequal(size(parameters),[1 4]));
assert(isequal(size(standard_deviation_in_z),[1 1]));
assert(isequal(size(z_fit),[Npoints 1]));
assert(isequal(size(unit_vector),[1 3]));
assert(isequal(size(standard_deviation_in_plane_orthogonals),[1 1]));

% Check variable values
assert(isequal(round(parameters,4),[0 0 1 0]));
assert(isequal(round(standard_deviation_in_z,4),0));
assert(isequal(round(z_fit,4),zeros(Npoints,1)));
assert(isequal(round(unit_vector,4),[0 0 1]));
assert(isequal(round(standard_deviation_in_plane_orthogonals,4),0));

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));

%% Compare speeds of pre-calculation versus post-calculation versus a fast variant
fig_num = 9903;
figure(fig_num);
close(fig_num);
rng(1823);

points = [0 0 0; 1 0 0; 0 1 0];

% Perform the calculation in slow mode
fig_num = [];
REPS = 100; minTimeSlow = Inf; 
tic;
for i=1:REPS
    tstart = tic;
    [parameters, standard_deviation_in_z, z_fit, unit_vector, base_point, standard_deviation_in_plane_orthogonals] = fcn_geometry_findPlaneNormal(points,fig_num);

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
    [parameters, standard_deviation_in_z, z_fit, unit_vector, base_point, standard_deviation_in_plane_orthogonals] = fcn_geometry_findPlaneNormal(points,fig_num);
    telapsed = toc(tstart);
    minTimeFast = min(telapsed,minTimeFast);
end
averageTimeFast = toc/REPS;

fprintf(1,'\n\nComparison of fast and slow modes of fcn_geometry_findPlaneNormal:\n');
fprintf(1,'N repetitions: %.0d\n',REPS);
fprintf(1,'Slow mode average speed per call (seconds): %.5f\n',averageTimeSlow);
fprintf(1,'Slow mode fastest speed over all calls (seconds): %.5f\n',minTimeSlow);
fprintf(1,'Fast mode average speed per call (seconds): %.5f\n',averageTimeFast);
fprintf(1,'Fast mode fastest speed over all calls (seconds): %.5f\n',minTimeFast);
fprintf(1,'Average ratio of fast mode to slow mode (unitless): %.3f\n',averageTimeSlow/averageTimeFast);
fprintf(1,'Fastest ratio of fast mode to slow mode (unitless): %.3f\n',minTimeSlow/minTimeFast);
%% Fail conditions
if 1==0
    %% FAIL 1: points not long enough
    points = [2 3];
    [root_point, unit_vector] = fcn_geometry_findPlaneNormal(points,fig_num);
    fprintf(1,'\n\nRoot point is: %.2f %.2f, Unit vector is: %.2f %.2f\n',root_point(1,1),root_point(1,2),unit_vector(1,1),unit_vector(1,2));
end