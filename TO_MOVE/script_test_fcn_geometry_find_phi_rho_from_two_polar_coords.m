% script_test_fcn_geometry_find_phi_rho_from_two_polar_coords
start_point_cart = [1 1];
end_point_cart   = [4 3];

[theta1,r1] = cart2pol(start_point_cart(:,1),start_point_cart(:,2));
[theta2,r2]   = cart2pol(end_point_cart(:,1),end_point_cart(:,2));

% Calculate the line
[phi,rho] = fcn_geometry_find_phi_rho_from_two_polar_coords(r1,r2,theta1,theta2);

%% Plot the results for points on the line
figure(1);
clf;
hold on;
grid minor;

plot([start_point_cart(1,1) end_point_cart(1,1)],[start_point_cart(1,2) end_point_cart(1,2)]);
angles = (min(theta1,theta2):0.01:max(theta1,theta2))';
radii_on_segment = rho./cos(angles - phi);
[x_points,y_points] = pol2cart(angles,radii_on_segment);
plot(x_points,y_points,'rx');

%% Plot the results for points around the line
figure(2);
clf;
hold on;
grid minor;

plot([start_point_cart(1,1) end_point_cart(1,1)],[start_point_cart(1,2) end_point_cart(1,2)]);

random_points = 5*rand(10,2);
plot(random_points(:,1),random_points(:,2),'k.');


[random_thetas,~] = cart2pol(random_points(:,1),random_points(:,2));

radii_on_segment = rho./cos(random_thetas - phi);
[x_points,y_points] = pol2cart(random_thetas,radii_on_segment);
plot(x_points,y_points,'rx');
