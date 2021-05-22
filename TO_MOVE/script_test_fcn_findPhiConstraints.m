%script_test_fcn_phi_constrainter
fig_num = 12121;

% Simple test - should result in 90 to 0
p_apex = [0 0];
vertex_1 = [1  0];
vertex_2 = [0 1];
[phi_start,change] = fcn_findPhiConstraints(p_apex,vertex_1,vertex_2,fig_num);
fprintf(1,'Results of phi constrainter: \n');
mod([phi_start change],2*pi)*180/pi


% Simple test - different cross product direction - should result in 90 to 0
p_apex = [0 0];
vertex_2 = [1  0];
vertex_1 = [0 1];
[phi_start,change] = fcn_findPhiConstraints(p_apex,vertex_1,vertex_2,fig_num);
fprintf(1,'Results of phi constrainter: \n');
mod([phi_start change],2*pi)*180/pi


% Simple test - should result in 0 to 90
p_apex = [0 0];
vertex_1 = [-1  0];
vertex_2 = [0 -1];
[phi_start,change] = fcn_findPhiConstraints(p_apex,vertex_1,vertex_2,fig_num);
fprintf(1,'Results of phi constrainter: \n');
mod([phi_start change],2*pi)*180/pi

% Simple test - should result in 90 to 45
p_apex = [0 0];
vertex_1 = [1  0];
vertex_2 = [-1 1];
[phi_start,change] = fcn_findPhiConstraints(p_apex,vertex_1,vertex_2,fig_num);
fprintf(1,'Results of phi constrainter: \n');
mod([phi_start change],2*pi)*180/pi

% Simple test to check the vectorization - should result in 0 to 90, and
% -45 to 45, 180 to 270
p_apex =   [0 0;  2                      0; 0 -1; -1 0];
vertex_1 = [1  0; 2+1/(2^0.5)  -1/(2^0.5) ; -1 -1; -2 0];
vertex_2 = [0  1; 2+1/(2^0.5)   1/(2^0.5) ; 0  -2; -1 1];

[phi_start,change] = fcn_findPhiConstraints(p_apex,vertex_1,vertex_2,fig_num);
fprintf(1,'Results of phi constrainter: \n');
mod([phi_start change],2*pi)*180/pi


angles = 0:0.01:2*pi;
for i_angle = 1:length(angles)
    angle = angles(i_angle);
    R = [cos(angle) -sin(angle); sin(angle) cos(angle)];
    
    p_apex2 = p_apex*R;
    vertex_12   = vertex_1*R;
    vertex_22   = vertex_2*R;
    
    [phi_start,change] = fcn_findPhiConstraints(p_apex2,vertex_12,vertex_22,fig_num);
    drawnow;
    %pause;
end

