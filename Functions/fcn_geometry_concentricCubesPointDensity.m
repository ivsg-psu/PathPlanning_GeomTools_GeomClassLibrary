function [points] = fcn_geometry_concentricCubesPointDensity(...
    exterior_size,interior_size,exterior_density,interior_density,varargin)

% Given the number of points and the range of the points, generates two
% concentric squares with the external having a length of (ext_length) and
% the interior square having a side length of half the ext_length rounded.
% Additional inputs can be given to add noise and a diagonal division line.
%
% FORMAT:
%
% [points] = ...
% fcn_geometry_concentricSquaresPointDensity(exterior_size,interior_size,
% exterior_density, interior_density,(noise),(diagonal_flag),(fig_num))
%
% INPUTS:
%   
%       exterior_size: Length of external side of the cube
%       interior_size: Length of internal side of the cube
%       exterior_density: Density of exterior area of the cube
%       interior_density: Density of interior area of the cube

%       (OPTIONAL INPUTS)
%
%       noise: Noise to give the figure
%       diagonal_flag: 1 or 0 input to have a diagonal half have noise
%       fig_num: Assigns a custom number to the figure
%
% OUTPUTS:
%       
%       points: An array of the X,Y,and Z positions of each point
%
% DEPENDENCIES:
%
%   None
%
% EXAMPLES: 
%
%       See the script:
%       script_test_fcn_geometry_concentricCubesPointDensity for a full
%       test suite.
%
% This function was written on 2024_6_15 by Aleksandr Goncharov
% Questions or comments? opg5041@psu.edu or 267-304-8354
%
% REVISIONS:
%
% 2024_6_20 - Aleksandr Goncharov
% Began rewriting code to become a cube, where density can be configured of
% internal and external areas.
%

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
    MATLABFLAG_LOADWZ_FLAG_CHECK_INPUTS = getenv("MATLABFLAG_LOADWZ_FLAG_CHECK_INPUTS");
    MATLABFLAG_LOADWZ_FLAG_DO_DEBUG = getenv("MATLABFLAG_LOADWZ_FLAG_DO_DEBUG");
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

% Calculate the number of points needed based on the size and density
extPoints = round(exterior_density * exterior_size^3);
intPoints = round(interior_density * interior_size^3);

% Calculate the number of points per axis
n_ext = round(nthroot(extPoints, 3));
n_int = round(nthroot(intPoints, 3));

% Recalculate spacing needed between points
spacing_ext = exterior_size / n_ext;
spacing_int = interior_size / n_int;

% Calculate side lengths from center
ext_side = exterior_size / 2;
int_side = interior_size / 2;

% Generate a grid of points spaced equidistant from each other
[x_e, y_e, z_e] = meshgrid(linspace(-ext_side, ext_side, n_ext), ...
                           linspace(-ext_side, ext_side, n_ext), ...
                           linspace(-ext_side, ext_side, n_ext));
[x_i, y_i, z_i] = meshgrid(linspace(-int_side, int_side, n_int), ...
                           linspace(-int_side, int_side, n_int), ...
                           linspace(-int_side, int_side, n_int));

% Flatten the matrix into a single column
x_ext = x_e(:);
y_ext = y_e(:);
z_ext = z_e(:);

x_int = x_i(:);
y_int = y_i(:);
z_int = z_i(:);

%remove points inside the internal zone
interior_mask=abs(x_ext) <= int_side & abs(y_ext) <= int_side & abs(z_ext) <= int_side;

x_ext=x_ext(~interior_mask);
y_ext=y_ext(~interior_mask);
z_ext=z_ext(~interior_mask);

%define placeholder displacement
x_ext_disp=zeros(size(x_ext));
y_ext_disp=zeros(size(x_ext));
z_ext_disp=zeros(size(x_ext));
x_int_disp=zeros(size(x_int));
y_int_disp=zeros(size(x_int));
z_int_disp=zeros(size(x_int));
  
%diagonal flag

%If diagonal 
if flag_create_diagonal == 1

    if noise == 0
    x_ext=x_ext+x_ext_disp;
    y_ext=y_ext+y_ext_disp;
    z_ext=z_ext+z_ext_disp;

    x_int=x_int+x_int_disp;
    y_int=y_int+y_int_disp;
    z_int=z_int+z_int_disp;
    end


    if noise ~= 0
        diag_ext_area = y_ext <= x_ext;
        diag_int_area = y_int <= x_int;

        x_ext_disp = (2*rand(size(x_ext(diag_ext_area),1),1) - 1) * spacing_ext * noise;
        y_ext_disp = (2*rand(size(y_ext(diag_ext_area),1),1) - 1) * spacing_ext * noise;
        z_ext_disp = (2*rand(size(z_ext(diag_ext_area),1),1) - 1) * spacing_ext * noise;

        x_int_disp = (2*rand(size(x_int(diag_int_area),1),1) - 1) * spacing_int * noise;
        y_int_disp = (2*rand(size(y_int(diag_int_area),1),1) - 1) * spacing_int * noise;
        z_int_disp = (2*rand(size(z_int(diag_int_area),1),1) - 1) * spacing_int * noise;


        x_ext(diag_ext_area)=x_ext(diag_ext_area)+x_ext_disp;
        y_ext(diag_ext_area)=y_ext(diag_ext_area)+y_ext_disp;
        z_ext(diag_ext_area)=z_ext(diag_ext_area)+z_ext_disp;

        x_int(diag_int_area)=x_int(diag_int_area)+x_int_disp;
        y_int(diag_int_area)=y_int(diag_int_area)+y_int_disp;
        z_int(diag_int_area)=z_int(diag_int_area)+z_int_disp;
    end

end

%If no diagonal flag is present
if flag_create_diagonal == 0
    if noise ~= 0
        x_ext_disp = (2*rand(size(x_ext,1),1) - 1) * spacing_ext * noise;
        y_ext_disp = (2*rand(size(y_ext,1),1) - 1) * spacing_ext * noise;
        z_ext_disp = (2*rand(size(z_ext,1),1) - 1) * spacing_ext * noise;

        x_int_disp = (2*rand(size(x_int,1),1) - 1) * spacing_int * noise;
        y_int_disp = (2*rand(size(y_int,1),1) - 1) * spacing_int * noise;
        z_int_disp = (2*rand(size(z_int,1),1) - 1) * spacing_int * noise;
    end


    x_ext=x_ext+x_ext_disp;
    y_ext=y_ext+y_ext_disp;
    z_ext=z_ext+z_ext_disp;

    x_int=x_int+x_int_disp;
    y_int=y_int+y_int_disp;
    z_int=z_int+z_int_disp;

end

%combining lists
X=[x_ext;x_int];
Y=[y_ext;y_int];
Z=[z_ext;z_int];

%making it as a point
points = [X Y Z];

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
    scatter3(x_ext, y_ext, z_ext, 'filled','c');
    scatter3(x_int, y_int, z_int, 'filled','r');
    title('Concentric Cubes Point Density Plot');
    axis equal;
    view(3)
    grid on;
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

end