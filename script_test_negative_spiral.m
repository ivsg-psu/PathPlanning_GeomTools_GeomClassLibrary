%%
close all

% Try positive curvature

arc1_parameters = [0    0.80000    0.8000   -3.5046   -1.9338         0    1.0000];
arc2_parameters = [0.2449    0.8000    0.8000   -1.9338   -0.3630         0    1.0000];
h0 = 0*pi/180;

[~, x0, y0, K0, Kf] = fcn_INTERNAL_calcArcParameters(h0, arc1_radius, arc2_radius);

spiralLength = 1;
s  = (0:0.01:1)'*spiralLength;

figure(1234);
clf;
hold on;
grid on;
axis equal;

% fcn_geometry_plotGeometry('arc',arc1_parameters);
% fcn_geometry_plotGeometry('arc',arc2_parameters);
[x_spiral,y_spiral] = fcn_geometry_extractXYfromSTSpiral(s,spiralLength,h0,x0,y0,K0,Kf,(1234));



% Find the center of the circle tangent at the end of the spiral
% Find the unit vector (need to do this analytically!)
s_tangent = [0.99999999 1]'*spiralLength;
[x_tangent,y_tangent] = fcn_geometry_extractXYfromSTSpiral(s_tangent,spiralLength,h0,x0,y0,K0,Kf);
unit_tangent = fcn_geometry_calcUnitVector([diff(x_tangent) diff(y_tangent)]);
approximate_end_angle = atan2(unit_tangent(2),unit_tangent(1))

analytical_end_angle   = h0 + (Kf-K0)*spiralLength/2 + K0*spiralLength
unit_tangent           = [cos(analytical_end_angle) sin(analytical_end_angle)];

unit_orthogonal = unit_tangent*[0 1; -1 0];
calculated_circle_center = arc2_radius*unit_orthogonal + [x_spiral(end) y_spiral(end)];



% Plot the circle's centers
plot(arc1_parameters(1,1),arc1_parameters(1,2),'b+');
plot(calculated_circle_center(:,1),calculated_circle_center(:,2),'r+');

% Plot the circles
fcn_geometry_plotCircle(arc1_parameters(1,1:2), arc1_parameters(1,3),'b-',(1234));
fcn_geometry_plotCircle(calculated_circle_center, arc2_radius,'r-',(1234));


%%
% Try negative curvature
arc1_parameters = [0    3.0000    3.0000   -2.9379   -1.3671         0    1.0000];
arc2_parameters = [0.6782   -0.8000    0.8000    1.7745    0.2037         0   -1.0000];


h0 = 0*pi/180;

[~, x0, y0, K0, Kf] = fcn_INTERNAL_calcArcParameters(h0, arc1_radius, -arc2_radius);

spiralLength = 1;
s  = (0:0.01:1)'*spiralLength;

figure(1234);
clf;
hold on;
grid on;
axis equal;

fcn_geometry_plotGeometry('arc',arc1_parameters);
% fcn_geometry_plotGeometry('arc',arc2_parameters);
[x_spiral,y_spiral] = fcn_geometry_extractXYfromSTSpiral(s,spiralLength,h0,x0,y0,K0,Kf,(1234));



% Find the center of the circle tangent at the end of the spiral
% Find the unit vector (need to do this analytically!)
s_tangent = [0.99999999 1]'*spiralLength;
[x_tangent,y_tangent] = fcn_geometry_extractXYfromSTSpiral(s_tangent,spiralLength,h0,x0,y0,K0,Kf);
unit_tangent = fcn_geometry_calcUnitVector([diff(x_tangent) diff(y_tangent)]);
approximate_end_angle = atan2(unit_tangent(2),unit_tangent(1))

analytical_end_angle   = h0 + (Kf-K0)*spiralLength/2 + K0*spiralLength
unit_tangent           = [cos(analytical_end_angle) sin(analytical_end_angle)];

unit_orthogonal = unit_tangent*[0 1; -1 0];
calculated_circle_center = -arc2_radius*unit_orthogonal + [x_spiral(end) y_spiral(end)];



% Plot the circle's centers
plot(arc1_parameters(1,1),arc1_parameters(1,2),'b+');
plot(calculated_circle_center(:,1),calculated_circle_center(:,2),'r+');

% Plot the circles
fcn_geometry_plotCircle(arc1_parameters(1,1:2), arc1_parameters(1,3),'b-',(1234));
fcn_geometry_plotCircle(calculated_circle_center, arc2_radius,'r-',(1234));


function [arc1_center_xy, x0, y0, K0, Kf] = fcn_INTERNAL_calcArcParameters(h0, arc1_radius, arc2_radius)
angle_in_arc1 = h0 - 90*pi/180;
arc1_center_xy = [0 arc1_radius];
arc1_xy_point  = arc1_center_xy + arc1_radius*[cos(angle_in_arc1) sin(angle_in_arc1)];
x0 = arc1_xy_point(1,1);
y0 = arc1_xy_point(1,2);
K0 = 1/arc1_radius;
Kf = 1/arc2_radius;
end
