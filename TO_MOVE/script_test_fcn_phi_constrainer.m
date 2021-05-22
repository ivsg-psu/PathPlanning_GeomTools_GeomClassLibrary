%script_test_fcn_phi_constrainter
fig_num = 12121;

% Simple test - should result in 0 to 90
p_apex = [0 0];
vertex_1 = [1  0];
vertex_2 = [0 1];
[phi_min,phi_max] = fcn_phi_constrainer(p_apex,vertex_1,vertex_2,fig_num);
fprintf(1,'Results of phi constrainter: \n');
mod([phi_min phi_max],2*pi)*180/pi

% Simple test - should result in 0 to 90
p_apex = [0 0];
vertex_1 = [-1  0];
vertex_2 = [0 -1];
[phi_min,phi_max] = fcn_phi_constrainer(p_apex,vertex_1,vertex_2,fig_num);
fprintf(1,'Results of phi constrainter: \n');
mod([phi_min phi_max],2*pi)*180/pi



% Simple test to check the vectorization - should result in 0 to 90, and
% -45 to 45, 180 to 270
p_apex =   [0 0;  2                      0; 0 -1];
vertex_1 = [1  0; 2+1/(2^0.5)  -1/(2^0.5) ; -1 -1];
vertex_2 = [0  1; 2+1/(2^0.5)   1/(2^0.5) ; 0  -2];

[phi_min,phi_max] = fcn_phi_constrainer(p_apex,vertex_1,vertex_2,fig_num);
fprintf(1,'Results of phi constrainter: \n');
mod([phi_min phi_max],2*pi)*180/pi


