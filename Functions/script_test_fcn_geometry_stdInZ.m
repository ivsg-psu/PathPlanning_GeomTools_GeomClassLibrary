%script test of: fcn_geometry_stdInZ
%Written by Aleksandr Goncharov

%% Basic Example
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   ____            _        ______                           _
%  |  _ \          (_)      |  ____|                         | |
%  | |_) | __ _ ___ _  ___  | |__  __  ____ _ _ __ ___  _ __ | | ___
%  |  _ < / _` / __| |/ __| |  __| \ \/ / _` | '_ ` _ \| '_ \| |/ _ \
%  | |_) | (_| \__ \ | (__  | |____ >  < (_| | | | | | | |_) | |  __/
%  |____/ \__,_|___/_|\___| |______/_/\_\__,_|_| |_| |_| .__/|_|\___|
%                                                      | |
%                                                      |_|
% See: https://patorjk.com/software/taag/#p=display&f=Big&t=Basic%20Example
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%§

% RUN PREVIOUS STEPS OF surfaceAnalysis, hard to emulate data otherwise

[mean_std_in_z_driven_path,mean_std_in_z_not_driven_path,max_std_in_z_not_driven_path]...
    = fcn_geometry_stdInZ(LiDAR_allPoints,updated_original_mapped_grids,gridIndices,gridIndices_cell_array,...
    gridCenters_updated_original_mapped_grids, updated_current_mapped_grids, current_grid_numbers_of_driven_path, ...
    gridCenters_driven_path,grid_AABBs, grid_size,1,2,3);

assert(mean_std_in_z_not_driven_path > mean_std_in_z_driven_path);
