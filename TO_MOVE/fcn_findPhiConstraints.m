function [phi_start,change_in_phi] = fcn_findPhiConstraints(p_apex,vertex_1,vertex_2,varargin)
%......................................................................
% [phi_min,phi_max] = fcn_findPhiConstraints(p_apex,vertex_1,vertex_2)
% This routine finds the minimum and maximum angle between the apex and
% centroid points given the path geometry
% inputs
% p_apex: [N x 2] apex points in an x,y format
% vertex_1: [N x 2] first adjacent vertex to be encountered by path
% vertex_2: [N x 2] second adjacent vertex to be encountered by path
% outputs
% phi_min: given in radians, the lower limit for the range of angles
% centroid points can be wrt the apex points
% phi_max: given in radians, the upper limit for the range of angles
% centroid points can be wrt the apex points
%  2020_04_06  -- written by Veronica Gruning (vag5076@psu.edu)
%  2020_05_02  - edits by Sean Brennan to add more comments, check results,
%  add debug option
%........................................................................

do_debug = 0;

%% check input arguments to see if we need to make a figure
if nargin < 3 || nargin > 4
    error('Incorrect number of input arguments. Expecting 3 or 4 arguments.')
end

if 4 == nargin % Did user give a figure number?
    fig_num = varargin{1};
    figure(fig_num);
    do_debug = 1; % To force plotting
else
    if do_debug  % Only create a new figure if debugging
        fig = figure(1); %#ok<UNRCH> % create new figure with next default index
        fig_num = get(fig,'Number');
    end
end

%% Code starts here
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

if do_debug
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
end
