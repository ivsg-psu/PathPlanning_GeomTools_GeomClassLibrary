function LLA_data = fcn_geometry_plotGridCenters(ENU_data,marker_size,RGB_triplet,varargin)
% This function takes in various boundary points and then generates a plot
% that shows where the boundary,driveable/non-driveable area, and unmapped area
% 
% FORMAT:  
% LLA_data =
% fcn_geometry_plotGridCentersBoundaryPoints(ENU_data,marker_size,RGB_triplet,(fig_num))
%
% INPUTS:
% ENU_data: Data points 
% marker_size: size of the points that are going to be plotted
% RGB_triplet: Color of the markers
% (OPTIONAL INPUT)
% legend_options: enable the plot to display a legend
% legend_name: name of the legend
% legend_position: position of the legend
% fig_num: figure number
%
% OUTPUTS:
% LLA coordinates of each point
%
% DEPENDENCIES:
% GPS CLASS
%
% EXAMPLES:
% See script: script_test_fcn_geometry_plotGridCenters
%
% Revision History
% 2024_07_10 - Aneesh Batchu 
% -- wrote the code originally
% 2024_07_10 - Aleksandr Goncharov
% -- Funtionlized this code
% 2024_07_11 - Aleksandr Goncharov
% -- Shortened and cleaned up the function to be more universal
% 2024_07-12 - Aneesh Batchu
% -- Changed the name of the function to "fcn_geometry_plotGridCenters"

flag_do_debug = 0; % Flag to plot the results for debugging
flag_max_speed = 0; % Flag to perform input checking

if flag_do_debug
    st = dbstack; %#ok<*UNRCH>
    fprintf(1,'STARTING function: %s, in file: %s\n',st(1).name,st(1).file);
    debug_fig_num = 34838; 
else
    debug_fig_num = []; 
end
%
%% check input arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____                   _
%  |_   _|                 | |
%    | |  _ __  _ __  _   _| |_ ___
%    | | | '_ \| '_ \| | | | __/ __|
%   _| |_| | | | |_) | |_| | |_\__ \
%  |_____|_| |_| .__/ \__,_|\__|___/
%              | |
%              |_|
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if flag_max_speed==0
    % Are there the right number of inputs?
    narginchk(3,7);
end

%does user want a legend
flag_create_legend=0;
if 4<=nargin
    temp=varargin{1};
    if temp==1
        flag_create_legend=1;
    end
end

%legend name
if 5<=nargin
    temp=varargin{2};
    if ~isempty(temp)
        legend_name=temp;
    end
end

%legend position
if 6<=nargin
    temp=varargin{3};
    if ~isempty(temp)
        legend_position=temp;
    end
end



% Does user want to specify fig_num?
fig_num = []; % Default is to have no figure
flag_do_plots = 0;
if (0==flag_max_speed) && (7<= nargin)
    temp = varargin{end};
    if ~isempty(temp)
        fig_num = temp;
        flag_do_plots = 1;
    end
end


%% Write main code for plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _
%  |  \/  |     (_)
%  | \  / | __ _ _ _ __
%  | |\/| |/ _` | | '_ \
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % True boundary points
% Trace_coordinates = [true_boundary_points,zeros(length(true_boundary_points),1)]; 

 % Define GPS object
reference_latitude = 40.86368573;
reference_longitude = -77.83592832;
reference_altitude = 344.189;
gps_object = GPS(reference_latitude,reference_longitude,reference_altitude); % Load the GPS class

% Use the class to convert LLA to ENU
LLA_data = gps_object.ENU2WGSLLA(ENU_data);


%% Any debugging?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____       _
%  |  __ \     | |
%  | |  | | ___| |__  _   _  __ _
%  | |  | |/ _ \ '_ \| | | |/ _` |
%  | |__| |  __/ |_) | |_| | (_| |
%  |_____/ \___|_.__/ \__,_|\__, |
%                            __/ |
%                           |___/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if flag_do_plots

% Plot the LLA boundary points
figure(fig_num);

% Plot
geoplot(LLA_data(:,1),LLA_data(:,2),'.','MarkerSize',marker_size,'Color',RGB_triplet,'DisplayName',legend_name);
hold on

%Adding Legend
if flag_create_legend == 1
    legend('Position',legend_position);
end

title('Boundary Points in LLA ')
geobasemap satellite
geotickformat -dd

end
end
