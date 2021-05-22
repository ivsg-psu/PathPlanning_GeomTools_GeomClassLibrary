
function [phi,rho] = fcn_geometry_find_phi_rho_from_two_polar_coords(r1,r2,theta1,theta2)
% Converts two polar points into the polar form of a line, giving phi and
% rho values for the line.
% See: http://www.nabla.hr/Z_MemoHU-015.htm
yterms = r1*cos(theta1)-r2*cos(theta2);
xterms = r2*sin(theta2)-r1*sin(theta1);
phi = atan2(yterms,xterms);
rho = r1*cos(theta1-phi);
end

