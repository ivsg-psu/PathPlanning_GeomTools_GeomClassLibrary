% fcn_geometry_plotLidarDeviation
% Write a function which takes x,y,z data,color_map,three figure numbers as the input
% inputs: cell_start, cell_end,color_map,fig_1,fig_2,fig_3)
% outputs: z, mean of z, angles, mean of angles

% Revision History
% Aneesh Batchu - 2024_07_08
% -- wrote the code originally
% Aleksandr Goncharov 2024_07_08
% -- functionalized the code


%Example:[z,mean_z,deviation_from_mean,angles,mean_angle,angle_difference] = fcn_geometry_plotLidarDeviation(50,59,'jet',111,222,333);

%%
function [z,mean_z,deviation_from_mean,angles,mean_angle,angle_difference] = fcn_geometry_plotLidarDeviation(cell_start,cell_end,color_map,fig_1,fig_2,fig_3)

load('LiDAR_ENU_Cell_Outer_Edge.mat','LiDAR_ENU_cell')
% LiDAR data
LiDAR_outer_edge = vertcat(LiDAR_ENU_cell{cell_start:cell_end});

x = LiDAR_outer_edge(:,1);
y = LiDAR_outer_edge(:,2);
z = LiDAR_outer_edge(:,3); 

% Calculate the deviation of z-values from the mean
mean_z = mean(z);
deviation_from_mean = z - mean_z;

figure(fig_1)
% Create the scatter plot
scatter3(x, y, z, 20, deviation_from_mean, 'filled');
colormap(color_map); 
colorbar; 
view(2); 
title('Deviation of z-value from the mean');
xlabel('x');
ylabel('y');

%%
% Calculate the angles relative to the x-axis
angles = atan2(y, x);

% Calculate the difference between each angle and the mean angle
mean_angle = mean(angles);
angle_difference = angles - mean_angle;

figure(fig_2)
% Create the scatter plot
scatter3(x, y, z, 20, angle_difference, 'filled');
colormap(color_map); 
colorbar; 
view(2); 
title('Angle difference from the mean angle');
xlabel('x');
ylabel('y');

%%

figure(fig_3)
colors = z;
  scatter3(x,y,z,20,colors);
colormap(color_map); % You can use 'viridis' if you have it, or other colormaps
        colorbar; % Add a color bar to show the color mapping
        view(2)
end