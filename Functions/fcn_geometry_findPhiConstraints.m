function [phi_start,change_in_phi] = ...
    fcn_geometry_findPhiConstraints(...
    p_apex,...
    vertex_1,...
    vertex_2,...
    varargin)
% fcn_geometry_findPhiConstraints
% This routine finds the minimum and maximum angle between the apex and
% centroid points given the path geometry
%
% FORMAT:
%
% function [phi_start,change_in_phi] = ...
%     fcn_geometry_findPhiConstraints(...
%     p_apex,...
%     vertex_1,...
%     vertex_2,...
%     (fig_num))
%
% INPUTS:
%
%      p_apex: [N x 2] list of apex points in an x,y format
%
%      vertex_1: [N x 2] list of first adjacent vertex to be encountered by path
%
%      vertex_2: [N x 2] list of second adjacent vertex to be encountered by path
%
%      (OPTIONAL INPUTS)
%
%      fig_num: a figure number to plot results.
%
% OUTPUTS:
%
%      phi_min: given in radians, the lower limit for the range of angles
%      that centroid points can be wrt the apex points 
%
%      phi_max: given in radians, the upper limit for the range of angles
%      that centroid points can be wrt the apex points
%
%
% DEPENDENCIES:
%
%      fcn_geometry_checkInputsToFunctions
%      fcn_geometry_plotCircle
%
% EXAMPLES:
%
% See the script: script_test_fcn_geometry_findPhiConstraints
% for a full test suite.
%
% This function was written on 2020_04_06 by V. Gruning
% Questions or comments? sbrennan@psu.edu

% Revision History:
%  2020_04_06  
%  -- written by Veronica Gruning (vag5076@psu.edu)
%  2020_05_02  
%  -- edits by Sean Brennan to add more comments, check results,
%  add debug option
%  2021_05_27  
%  -- edits for input checking, prep for geometry library


%% Debugging and Input checks
flag_check_inputs = 1; % Set equal to 1 to check the input arguments
flag_do_plot = 0;      % Set equal to 1 for plotting
flag_do_debug = 0;     % Set equal to 1 for debugging

if flag_do_debug
    st = dbstack; %#ok<*UNRCH>
    fprintf(1,'STARTING function: %s, in file: %s\n',st(1).name,st(1).file);
end

%% check input arguments?
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

%% check input arguments 
   
    
if flag_check_inputs    
    % Are there the right number of inputs?
    narginchk(3,4);
    
    % Check the p_apex input
    fcn_geometry_checkInputsToFunctions(...
        p_apex, '2column_of_numbers');
    
    % Use number of radii to calculate the number of centers
    N_apexes = length(p_apex(:,1));
    
    % Check the vertex_1 input
    fcn_geometry_checkInputsToFunctions(...
        vertex_1, '2column_of_numbers',[N_apexes N_apexes]);
    
    % Check the vertex_2 input
    fcn_geometry_checkInputsToFunctions(...
        vertex_2, '2column_of_numbers',[N_apexes N_apexes]);
    
end
    

% Does user want to show the plots?
if 4 == nargin
    fig_num = varargin{end};
    figure(fig_num);
    flag_do_plot = 1;    
else
    if flag_do_debug
        fig = figure;
        fig_num = fig.Number;
        flag_do_plot = 1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _
%  |  \/  |     (_)
%  | \  / | __ _ _ _ __
%  | |\/| |/ _` | | '_ \
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate vectors corresonding to apexes
p_apex_vertex_1 = vertex_1 - p_apex; % vector from apex to first vertex
p_apex_vertex_2 = vertex_2 - p_apex; % vector from apex to second vertex

% Calculate the start angle of vertex 1
angle_vertex_1 = atan2(p_apex_vertex_1(:,2), p_apex_vertex_1(:,1)); % angle of apex to first vertex

% Calculate the angle between the vectors
mag_1 = sum(p_apex_vertex_1.^2,2).^0.5;
mag_2 = sum(p_apex_vertex_2.^2,2).^0.5;
change = acos(sum(p_apex_vertex_1.*p_apex_vertex_2,2)./(mag_1.*mag_2));

% Calculate the cross product, to determine directions
cross_vertex_1_to_vertex_2 = ...
    cross(...
    [p_apex_vertex_1 zeros(length(p_apex_vertex_1(:,1)),1)],...
    [p_apex_vertex_2 zeros(length(p_apex_vertex_2(:,1)),1)],2);
% ... and keep just the z-part
is_positive_cross_vertex_1_to_vertex_2 = 2*(cross_vertex_1_to_vertex_2(:,3)>0)-1;

% Calculate the resulting change
% If cross product is positive, then segment 1 is sweeping clockwise to
% segment 2. The start angle will then be (angle1 + 90) degrees
% and the end angle will be (angle1 + change - 90 degrees). The difference
% is thus going to be (change - 180) degrees (e.g. in the negative direction).
%
% If the cross product is negative, then the angle starts at (angle1 - 90)
% and it ends at (angle1 - change + 90 degrees). The difference is then
% (180 - change) degrees, (e.g. in the positive direction).
phi_start = angle_vertex_1 + pi/2.*is_positive_cross_vertex_1_to_vertex_2;
change_in_phi = (change - pi).*is_positive_cross_vertex_1_to_vertex_2;


%% Plot results?
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

if flag_do_plot
    % Set up the figure
    figure(fig_num);
    clf;
    grid on;
    hold on;
    axis equal;

    % Plot points
    plot(p_apex(:,1),p_apex(:,2),'rx');
    text(p_apex(:,1),p_apex(:,2),'apex');
    
    plot(vertex_1(:,1),vertex_1(:,2),'go');
    text(vertex_1(:,1),vertex_1(:,2),'vertex_1');
    
    plot(vertex_2(:,1),vertex_2(:,2),'ro');
    text(vertex_2(:,1),vertex_2(:,2),'vertex_2');
    
    % Plot lines
    N_apexes = length(p_apex(:,1));
    
    % Lines for vertex_1
    vertex_lines_x = [p_apex(:,1) vertex_1(:,1) NaN*p_apex(:,1)];
    vertex_lines_y = [p_apex(:,2) vertex_1(:,2) NaN*p_apex(:,2)];
    vertex_lines_x = reshape(vertex_lines_x',N_apexes*3,1);
    vertex_lines_y = reshape(vertex_lines_y',N_apexes*3,1);   
    plot(vertex_lines_x,vertex_lines_y,'g-');
    
    % Lines for vertex_2
    vertex_lines_x = [p_apex(:,1) vertex_2(:,1) NaN*p_apex(:,1)];
    vertex_lines_y = [p_apex(:,2) vertex_2(:,2) NaN*p_apex(:,2)];
    vertex_lines_x = reshape(vertex_lines_x',N_apexes*3,1);
    vertex_lines_y = reshape(vertex_lines_y',N_apexes*3,1);   
    plot(vertex_lines_x,vertex_lines_y,'r-');
    
    % Show results of angle calculations
    middles_v1 = (p_apex+vertex_1)/2;
    middles_v2 = (p_apex+vertex_2)/2;
    middles_v1_to_v2 = (vertex_1+vertex_2)/2;
    middles_mid_to_apex = (middles_v1_to_v2+p_apex)/2;
    for i_apex = 1:N_apexes
        text(middles_v1(i_apex,1),middles_v1(i_apex,2),sprintf('%.f deg',angle_vertex_1(i_apex,1)*180/pi));
        text(middles_v2(i_apex,1),middles_v2(i_apex,2),sprintf('%.f deg',(angle_vertex_1(i_apex,1)+change(i_apex,1).*is_positive_cross_vertex_1_to_vertex_2(i_apex,1))*180/pi));
        text(middles_mid_to_apex(i_apex,1),middles_mid_to_apex(i_apex,2),sprintf('Searching from %.f deg with change of %.f deg',phi_start(i_apex,1)*180/pi,change_in_phi(i_apex,1)*180/pi));
    end
    
end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file); 
end


end % Ends the function
