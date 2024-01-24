function sphere_points = fcn_geometry_plotSphere(centers,radii,varargin)   
% fcn_geometry_plotSphere -  plots a sphere by creating a vector of angles
% spaced 0.01 radians apart, and plotting this as a line around the
% perimeter.
%
% FORMAT:
%
%     XYZ_points = fcn_geometry_plotSphere(...
%     centers,...
%     radii,...
%     (color_vector);
%     (fig_num))
%
% INPUTS:
%
%      centers: an [N x 2] vector in [x y] of the points of sphere centers
%
%      radii: a [N x 1] vector of the radii of the spheres (to avoid
%      calculation time)
%
%      (OPTIONAL INPUTS)
%
%      color_vector: A color vector, e.g. [1 0 0.23], that dictates the
%      sphere color. 
%
%      fig_num: a figure number to plot results.
%
% OUTPUTS:
%
%      sphere_points: the [x y z] coordinates of the sphere points. If N
%      centers and radii are given, with N>1, then sphere_points will be a
%      cell array of points.
%
% DEPENDENCIES:
%
%      fcn_DebugTools_checkInputsToFunctions
%
% EXAMPLES:
%
% See the script: script_test_fcn_geometry_plotSphere
% for a full test suite.
%
% This function was written on 2020_10_13 by S. Brennan
% Questions or comments? sbrennan@psu.edu

% Revision History:
% 2021-05-22
% -- new function from fcn_geometry_findAngleUsing3PointsOnCircle
% -- eliminates repo on fcn_plotCircles


%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==4 && isequal(varargin{end},-1))
    flag_do_debug = 0; % Flag to plot the results for debugging
    flag_check_inputs = 0; % Flag to perform input checking
    flag_max_speed = 1;
else
    % Check to see if we are externally setting debug mode to be "on"
    flag_do_debug = 0; % Flag to plot the results for debugging
    flag_check_inputs = 1; % Flag to perform input checking
    MATLABFLAG_GEOMETRY_FLAG_CHECK_INPUTS = getenv("MATLABFLAG_GEOMETRY_FLAG_CHECK_INPUTS");
    MATLABFLAG_GEOMETRY_FLAG_DO_DEBUG = getenv("MATLABFLAG_GEOMETRY_FLAG_DO_DEBUG");
    if ~isempty(MATLABFLAG_GEOMETRY_FLAG_CHECK_INPUTS) && ~isempty(MATLABFLAG_GEOMETRY_FLAG_DO_DEBUG)
        flag_do_debug = str2double(MATLABFLAG_GEOMETRY_FLAG_DO_DEBUG);
        flag_check_inputs  = str2double(MATLABFLAG_GEOMETRY_FLAG_CHECK_INPUTS);
    end
end


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
if (0==flag_max_speed)
    if flag_check_inputs
        % Are there the right number of inputs?
        narginchk(2,4);

        % Check the centers input
        fcn_DebugTools_checkInputsToFunctions(...
            centers, '3column_of_numbers');

        % Use number of radii to calculate the number of centers
        Nspheres = length(centers(:,1));

        % Check the radii input
        fcn_DebugTools_checkInputsToFunctions(...
            radii, '1column_of_numbers',Nspheres);

    end
end

% Does user want to specify the color_vector? 
color_vector = [];
if (3 <=nargin)
    temp = varargin{1};
    if ~isempty(temp)
        color_vector = temp;
    end
end

% Does user want to specify the figure? And make sure no plots are formed
% if on max_speed mode!
if (0==flag_max_speed)
    flag_do_plot = 1;

    if (4 <= nargin)
        temp = varargin{end};
        if ~isempty(temp)
            fig_num = temp;
        else
            flag_do_plot = 0;
        end
    else
        fig = gcf;
        fig_num = fig.Number;
    end
else
    flag_do_plot = 0;
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


% Get a default sphere construct
[X,Y,Z] = sphere;

% How many spheres are we calculating?
N_spheres = length(radii);

% Initialize cell array, if more than one sphere
if N_spheres~=1
    sphere_points{N_spheres} = [];
end

% Loop through spheres, calculating XYZ data for each
plotting_X{N_spheres}= [];
plotting_Y{N_spheres}= [];
plotting_Z{N_spheres}= [];

for ith_sphere = 1:N_spheres
    % Fix the size and position
    X2 = X * radii(ith_sphere) + centers(ith_sphere,1);
    Y2 = Y * radii(ith_sphere) + centers(ith_sphere,2);
    Z2 = Z * radii(ith_sphere) + centers(ith_sphere,3);


    % Save results
    plotting_X{ith_sphere}= X2;
    plotting_Y{ith_sphere}= Y2;
    plotting_Z{ith_sphere}= Z2;

    if nargout>0
        sizeX = size(X2);
        X3 = reshape(X2,sizeX(1)*sizeX(2),1);
        Y3 = reshape(Y2,sizeX(1)*sizeX(2),1);
        Z3 = reshape(Z2,sizeX(1)*sizeX(2),1);
        if N_spheres==1
            sphere_points = [X3 Y3 Z3];
        else
            sphere_points{ith_sphere} = [X3 Y3 Z3];
        end

    end
end

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
    % set up the figure and check if it has been used before. If not used
    % already, note that the axes will be re-scaled
    temp_h = figure(fig_num);
    flag_rescale_axis = 0;
    if isempty(get(temp_h,'Children'))
        flag_rescale_axis = 1;
        hold on;
        axis equal
        grid minor;
        view(3);
    else
        hold on;
        axis equal
    end

    for ith_sphere = 1:N_spheres

        X2 = plotting_X{ith_sphere};
        Y2 = plotting_Y{ith_sphere};
        Z2 = plotting_Z{ith_sphere};


        % Make plots
        if ~isempty(color_vector)
            N_levels = 10;
            increment = 1/N_levels;
            scalings = (increment:increment:1)';
            if length(color_vector(:,1))==length(radii)
                current_color_vector = color_vector(ith_sphere,:);
            else
                current_color_vector = color_vector ;
            end
            color_map = current_color_vector .* scalings;
            surf(X2,Y2,Z2,'EdgeColor',current_color_vector*increment,'FaceAlpha',0.1,'EdgeAlpha',0.1);

            colormap(color_map);
        else
            surf(X2,Y2,Z2,'EdgeColor',[0.4 0.4 0.4],'FaceAlpha',0.1,'EdgeAlpha',0.1);
        end

    end

    % Make axis slightly larger?
    if flag_rescale_axis
        temp = axis;
        %     temp = [min(points(:,1)) max(points(:,1)) min(points(:,2)) max(points(:,2))];
        axis_range_x = temp(2)-temp(1);
        axis_range_y = temp(4)-temp(3);
        percent_larger = 0.3;
        axis([temp(1)-percent_larger*axis_range_x, temp(2)+percent_larger*axis_range_x,  temp(3)-percent_larger*axis_range_y, temp(4)+percent_larger*axis_range_y]);

    end
end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file); 
end
