function [points] = fcn_geometry_concentricSquaresPointDensity_new(ext_length,int_length,ext_point_concentration,int_point_concentration,varargin)

% Given the number of points and the range of the points, generates two
% concentric squares with the external having a length of (ext_length) and
% the interior square having a side length of half the ext_length rounded.
% Additional inputs can be given to add noise and a diagonal division line.
%
% FORMAT:
%
% [X_coord,Y_coord,Z_coord] = ...
% fcn_geometry_concentricSquaresPointDensity(N_points, ext_length, ..
% (noise_lvl), (diagonal_flag),(fig_num))
%
% INPUTS:
%   
%       ext_length: Length of external side of the square
%       int_length: Length of internal side of the square
%       ext_point_concentration: Concentration of external points per unit^2
%       int_point_concentration: Concentration of internal points per unit^2
%       (OPTIONAL INPUTS)
%
%       noise: Noise to give the figure
%       diagonal_flag: 1 or 0 input to have a diagonal half have noise
%       fig_num: Assigns a custom number to the figure
%
% OUTPUTS:
%       
%       Points: Values of the X,Y,and Z positions of each point
%
% DEPENDENCIES:
%
%   None
%
% EXAMPLES: 
%
%       See the script:
%       script_test_fcn_geometry_concentricSquaresPointDensity_new for a full
%       test suite.
%
% This function was written on 2024_6_15 by Aleksandr Goncharov
% Questions or comments? opg5041@psu.edu or 267-304-8354
%
%Revision History:
% 2024_6_24 - A. Goncharov - Revised to have point concentration input instead of N_points

%% Debug and Max speed
% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.

flag_max_speed = 0;
if (nargin==7 && isequal(varargin{end},-1))
    flag_do_debug = 0; % Flag to plot the results for debugging
    flag_max_speed = 1;
else
    % Check to see if we are externally setting debug mode to be "on"
    flag_do_debug = 0; % Flag to plot the results for debugging
    MATLABFLAG_LOADWZ_FLAG_CHECK_INPUTS = getenv("MATLABFLAG_GEOMETRY_FLAG_CHECK_INPUTS");
    MATLABFLAG_LOADWZ_FLAG_DO_DEBUG = getenv("MATLABFLAG_GEOMETRY_FLAG_DO_DEBUG");
    if ~isempty(MATLABFLAG_LOADWZ_FLAG_CHECK_INPUTS) && ~isempty(MATLABFLAG_LOADWZ_FLAG_DO_DEBUG)
        flag_do_debug = str2double(MATLABFLAG_LOADWZ_FLAG_DO_DEBUG);
    end
end
if flag_do_debug
    st = dbstack; %#ok<*UNRCH>
    fprintf(1,'STARTING function: %s, in file: %s\n',st(1).name,st(1).file);
    debug_fig_num = 34838; %#ok<NASGU>
else
    debug_fig_num = []; %#ok<NASGU>
end




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

% flag_do_plots = 1;
% fig_num=1;

if flag_max_speed == 0
    % Are there the right number of inputs?
    narginchk(4,7);
end


if nargin==4
    noise=0;
   
end

%Did the user want noise?
if 5 <= nargin
    temp = varargin{1};
    noise=0;
    if ~isempty(temp)
        noise=temp;
    end
end

%Did the user want a diagonal?
flag_create_diagonal=0;
if 6 <= nargin
    temp = varargin{2};
    if temp == 1
        flag_create_diagonal = 1;
    end
end

% %Does user want specific fig_num?
% if 5 == nargin
%     temp = varargin{end};
%     if ~isempty(temp)
%         fig_num = temp;
%     end
% end

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

%Calculate number of points needed independently
ext_points= round(ext_point_concentration*ext_length^2);
int_points= round(int_point_concentration*int_length^2);

n_ext = round(sqrt(ext_points));
n_int = round(sqrt(int_points));

spacing_ext=ext_length/n_ext;
spacing_int=int_length/n_int;

ext_side=ext_length/2;
int_side=int_length/2;

%grid generation outputs X-Y coordinates of every point
[x_e,y_e]=meshgrid(linspace(-ext_side,ext_side,n_ext),linspace(-ext_side,ext_side,n_ext));
x_ext=x_e(:);
y_ext=y_e(:);

[x_i,y_i]=meshgrid(linspace(-int_side,int_side,n_int),linspace(-int_side,int_side,n_int));
x_int=x_i(:);
y_int=y_i(:);


%point remover, only keeps points outside the masked area
interior_mask= abs(x_ext) <= int_side & abs(y_ext) <= int_side;

x_ext=x_ext(~interior_mask);
y_ext=y_ext(~interior_mask);
z_ext=zeros(length(x_ext),1);
z_int=zeros(length(x_int),1);

if noise~=0

    x_ext_disp = (-spacing_ext + (spacing_ext*2)*rand(size(x_ext,1),1)) * noise;
    y_ext_disp = (-spacing_ext + (spacing_ext*2)*rand(size(y_ext,1),1)) * noise;
    z_ext_disp = (-spacing_ext + (spacing_ext*2)*rand(size(z_ext,1),1)) * noise;

    x_int_disp = (-spacing_int + (spacing_int*2)*rand(size(x_int,1),1)) * noise;
    y_int_disp = (-spacing_int + (spacing_int*2)*rand(size(y_int,1),1)) * noise;
    z_int_disp = (-spacing_int + (spacing_int*2)*rand(size(z_int,1),1)) * noise;

    x_ext=x_ext+x_ext_disp;
    y_ext=y_ext+y_ext_disp;
    z_ext=z_ext+z_ext_disp;

    x_int=x_int+x_int_disp;
    y_int=y_int+y_int_disp;
    z_int=z_int+z_int_disp;

end

%Diagonal Check
if flag_create_diagonal==1
    %Create diag line points
    diag_x=[-ext_side,ext_side];
    diag_y=[-ext_side,ext_side];

    %Create diagonal mask
    diag_area_ext= y_ext<=x_ext;
    diag_area_int= y_int<=x_int;

    if noise ~=0
    %Revert a diagonal area to an undisturbed area
    x_ext(diag_area_ext)=x_ext(diag_area_ext)-x_ext_disp(diag_area_ext);
    y_ext(diag_area_ext)=y_ext(diag_area_ext)-y_ext_disp(diag_area_ext);
    z_ext(diag_area_ext)=z_ext(diag_area_ext)-z_ext_disp(diag_area_ext);

    x_int(diag_area_int)=x_int(diag_area_int)-x_int_disp(diag_area_int);
    y_int(diag_area_int)=y_int(diag_area_int)-y_int_disp(diag_area_int);
    z_int(diag_area_int)=z_int(diag_area_int)-z_int_disp(diag_area_int);

    end
end

X=[x_ext;x_int];
Y=[y_ext;y_int];
Z=[z_ext;z_int];

points=[X Y Z];

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


%Plotting the Results

if flag_do_plots

    
    figure(fig_num);
    hold on
  

    if flag_create_diagonal == 1
        %plots the diagonal
        plot(diag_x, diag_y,'black' ,'LineWidth', 2);

    end


    scatter3(x_ext, y_ext, z_ext, 'filled','c');
    scatter3(x_int, y_int, z_int, 'filled','r');
    view(3)
    hold off

end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file);
end

%% Functions follow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   ______                _   _
%  |  ____|              | | (_)
%  | |__ _   _ _ __   ___| |_ _  ___  _ __  ___
%  |  __| | | | '_ \ / __| __| |/ _ \| '_ \/ __|
%  | |  | |_| | | | | (__| |_| | (_) | | | \__ \
%  |_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|___/
%
% See: https://patorjk.com/software/taag/#p=display&f=Big&t=Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%§

