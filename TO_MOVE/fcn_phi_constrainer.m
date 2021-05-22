function [phi_in,phi_out] = fcn_phi_constrainer(p_apex,vertex_1,vertex_2,varargin)
%......................................................................
% [phi_min,phi_max] = fcn_phi_constrainer(p_apex,vertex_1,vertex_2)
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

% Calculate the angles of the vectors relative to the x-axis
angle_vertex_1 = atan2(p_apex_vertex_1(:,2), p_apex_vertex_1(:,1)); % angle of apex to first vertex
angle_vertex_2 = atan2(p_apex_vertex_2(:,2), p_apex_vertex_2(:,1)); % angle of apex to second vertex

% Calculate the cross product, to determine directions
cross_vertex_1_to_vertex_2 = ...
    cross(...
    [p_apex_vertex_1 zeros(length(p_apex_vertex_1(:,1)),1)],...
    [p_apex_vertex_2 zeros(length(p_apex_vertex_2(:,1)),1)]);
% ... and keep just the z-part
cross_vertex_1_to_vertex_2 = cross_vertex_1_to_vertex_2(:,3);


phi_in = (angle_vertex_1 + pi/2).*cross_vertex_1_to_vertex_2 + (angle_vertex_1 - pi/2).*(1-cross_vertex_1_to_vertex_2);
phi_out = (angle_vertex_2 - pi/2).*cross_vertex_1_to_vertex_2 + (angle_vertex_2 + pi/2).*(1-cross_vertex_1_to_vertex_2);

% phi_min = (angle_vertex_2 - pi/2) .*(angle_vertex_1<angle_vertex_2) + (angle_vertex_1 - pi/2).*(1-angle_vertex_1<angle_vertex_2);
% phi_max = (angle_vertex_1 + pi/2).*(angle_vertex_1<angle_vertex_2) + (angle_vertex_2 + pi/2) .*(1-angle_vertex_1<angle_vertex_2);

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
    for i_apex = 1:N_apexes
        text(middles_v1(i_apex,1),middles_v1(i_apex,2),sprintf('%.f deg',angle_vertex_1(i_apex,1)*180/pi));
        text(middles_v2(i_apex,1),middles_v2(i_apex,2),sprintf('%.f deg',angle_vertex_2(i_apex,1)*180/pi));
        text(middles_v1_to_v2(i_apex,1),middles_v1_to_v2(i_apex,2),sprintf('Searching from %.f deg to %.f deg',phi_in(i_apex,1)*180/pi,phi_out(i_apex,1)*180/pi));
    end
    
    
end
end
