function [points] = fcn_geometry_concentricSquaresPointDensity(N_points,ext_length,int_length,varargin)

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
%       N_Points: Number of points to plot
%       ext_length: Length of external side of the square
%       inr_length: Length of internal side of the square
%       (OPTIONAL INPUTS)
%
%       noise_lvl: Noise to give the figure
%       diagonal_flag: 1 or 0 input to have a diagonal half have noise
%       fig_num: Assigns a custom number to the figure
%
% OUTPUTS:
%       
%       X,Y,Z: Values of the X,Y,and Z positions of each point
%
% DEPENDENCIES:
%
%   None
%
% EXAMPLES: 
%
%       See the script:
%       script_test_fcn_GeomTools_concentricSquaresPointDensity for a full
%       test suite.
%
% This function was written on 2024_6_15 by Aleksandr Goncharov
% Questions or comments? opg5041@psu.edu or 267-304-8354

%% Debug and Max speed
% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.

flag_max_speed = 0;
if (nargin==6 && isequal(varargin{end},-1))
    flag_do_debug = 0; % Flag to plot the results for debugging
    flag_check_inputs = 0; % Flag to perform input checking
    flag_max_speed = 1;
else
    % Check to see if we are externally setting debug mode to be "on"
    flag_do_debug = 0; % Flag to plot the results for debugging
    flag_check_inputs = 1; % Flag to perform input checking
    MATLABFLAG_LOADWZ_FLAG_CHECK_INPUTS = getenv("MATLABFLAG_LOADWZ_FLAG_CHECK_INPUTS");
    MATLABFLAG_LOADWZ_FLAG_DO_DEBUG = getenv("MATLABFLAG_LOADWZ_FLAG_DO_DEBUG");
    if ~isempty(MATLABFLAG_LOADWZ_FLAG_CHECK_INPUTS) && ~isempty(MATLABFLAG_LOADWZ_FLAG_DO_DEBUG)
        flag_do_debug = str2double(MATLABFLAG_LOADWZ_FLAG_DO_DEBUG);
        flag_check_inputs  = str2double(MATLABFLAG_LOADWZ_FLAG_CHECK_INPUTS);
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
    narginchk(3,6);
end


if nargin==3
    noise=0;
   
end

%Did the user want noise?
flag_create_noise = 0;
if 4 <= nargin
    temp = varargin{1};
    noise=0;
    if ~isempty(temp)
        flag_create_noise = 1;
        noise=temp;
    end
end

%Did the user want a diagonal?
flag_create_diagonal=0;
if 5 <= nargin
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
if (0==flag_max_speed) && (6<= nargin)
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

%Height without noise
Z=0;

%Create The Squares

A=ext_length/2;  %Exterior Length from center
B=int_length/2;   %Interior Square from Center (Half of exterior square)

%Generate Points Inside density ~~3 times higher
pointsInner=round(N_points*3/4);
pointsOuter=N_points-pointsInner;

%Points outside the interior square
X_out = (rand(pointsOuter, 1) - 0.5) * ext_length;
Y_out = (rand(pointsOuter, 1) - 0.5) * ext_length;

%Points inside the interior square
X_in = (rand(pointsInner, 1) - 0.5) * int_length;
Y_in = (rand(pointsInner, 1) - 0.5) * int_length;

X=[X_out;X_in];
Y=[Y_out;Y_in];

%Creates points without a diagonal 
if flag_create_diagonal == 0

    Z_in=zeros(pointsInner,1);
    Z_out=zeros(pointsOuter,1);

    if noise~=0
        Z_in = Z + noise * randn(pointsInner, 1);
        Z_out = Z + noise * randn(pointsOuter, 1);
        
    end
    
    Z=[Z_out;Z_in];

end

%If Diagonal Flag is checked 
if flag_create_diagonal == 1

    %Points to create a diagonal line across
    diag_x = [-A, A];
    diag_y = [-A, A];
        
    Z_in=zeros(pointsInner,1);
    Z_out=zeros(pointsOuter,1);
    
    if noise ~=0
    %Filter to determine which points are in a diagonal sector
    noisy_area_out = Y_out <= X_out;
    noisy_area_in  = Y_in  <= X_in;
    %Differentiating noisy and non noisy points
    %Noisy
    noisy_points_X_out=X_out(noisy_area_out);
    noisy_points_Y_out=Y_out(noisy_area_out);
    noisy_points_X_in=X_in(noisy_area_in);
    noisy_points_Y_in=Y_in(noisy_area_in);
    %Non-noisy
    undisturbed_X_out = X_out(~noisy_area_out);
    undisturbed_Y_out = Y_out(~noisy_area_out);
    undisturbed_X_in = X_in(~noisy_area_in);
    undisturbed_Y_in = Y_in(~noisy_area_in);
    %number of points inside and outside that are noisy/undisturbed
    num_Z_noisy_out=zeros(size(noisy_points_X_out));
    num_Z_noisy_in=zeros(size(noisy_points_X_in));
    num_Z_undisturbed_out=zeros(size(undisturbed_X_out));
    num_Z_undisturbed_in=zeros(size(undisturbed_X_in));
    %apply noise
    noisy_Z_out=num_Z_noisy_out + noise * randn(sum(noisy_area_out),1);
    noisy_Z_in=num_Z_noisy_in + noise * randn(sum(noisy_area_in),1);
    % Combine noisy and non-noisy points
    Z_out(noisy_area_out) = noisy_Z_out;
    Z_in(noisy_area_in) = noisy_Z_in;
    
    
   
    end

    Z=[Z_out;Z_in];
end 
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
    %plot3(X,Y,Z)


    P=fill3([-A -A A A -A],[-A A A -A -A],[0 0 0 0 0],'r',[-B -B B B -B],[-B B B -B -B],[0 0 0 0 0],'b');
    P(1).FaceAlpha=0.05;
    P(2).FaceAlpha=0.05;
    hold on

    if flag_create_diagonal == 1

        %plots the diagonal
        plot(diag_x, diag_y,'black' ,'LineWidth', 2);

        if noise ~= 0

            scatter3(X_in, Y_in, Z_in, 'filled', 'red');
            scatter3(X_out, Y_out, Z_out, 'filled', 'blue');
        end

    end


    if flag_create_diagonal == 0 | noise==0
        scatter3(X_in,Y_in,Z_in,'filled','red')
        scatter3(X_out,Y_out,Z_out,'filled','blue')
        hold off
    end


    view(3)

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

