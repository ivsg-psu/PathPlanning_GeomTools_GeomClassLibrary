%script_fcn_geometry_plotGridCenters
%script for fcn_geometry_plotGridCenters
%
%
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

ENU_data =[418.5194  230.3028 -344.2069;
    419.1155  229.8859 -344.2069;
    419.8043  229.5207 -344.2069;
    420.3828  229.0694 -344.2070;
    421.0570  228.7612 -344.2070;
    421.5912  228.2227 -344.2070;
    422.2313  227.8078 -344.2070;
    422.9040  227.4922 -344.2071;
    423.4936  227.0386 -344.2071;
    424.0961  226.5475 -344.2071;
    424.7295  226.1757 -344.2071];

ENU_data2=ENU_data+1;
ENU_data3=ENU_data-1;

marker_size=30;
RGB_triplet=[0.8 0.8 0.8];
RGB_triplet2=[1 1 0];
RGB_triplet3=[0 1 1];
fig_num=1111;
legend_option=1;
legend_name='Legend Test';
legend_name2='Legend Test 2';
legend_name3='Legend Test 3';

%fcn_geometry_plotGridCentersBoundaryPoints(ENU_data,marker_size,RGB_triplet,(legend_option),(legend_name),(fig_num))
LLA_data = fcn_geometry_plotGridCenters(ENU_data,marker_size,RGB_triplet,legend_option,legend_name,fig_num);
LLA_data2 = fcn_geometry_plotGridCenters(ENU_data2,marker_size,RGB_triplet2,legend_option,legend_name2,fig_num);
LLA_data3 = fcn_geometry_plotGridCenters(ENU_data3,marker_size,RGB_triplet3,legend_option,legend_name3,fig_num);

%%