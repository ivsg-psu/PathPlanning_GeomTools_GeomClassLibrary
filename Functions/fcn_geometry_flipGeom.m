function geomParameters_flipped = fcn_geometry_flipGeom(geomType,  geomParameters, varargin) 
%% fcn_geometry_flipGeom
% This function "flips" a geometry given its type and parameter set such
% that the start/end of the input geometry becomes the end/start of the
% output geometry. Spatially, the geometry does not move.
% 
% FORMAT: 
%
% geomParameters_flipped = fcn_geometry_flipGeom(geomType,  geomParameters, (fig_num)) 
% 
% INPUTS:
%
%      geomType: The string input of the geometry. Can be one of the
%      following strings:
%   
%           'circle', 'arc', 'line', 'line segment' or 'segment'
%           (Pending: 'spiral')
%
%      geomParameters: The parameters of the geometry as a
%      [1xN] vector.  See fcn_geometry_fillEmptyDomainStructure for details
%      on how to fill the vector for different geometry types.
%
% (OPTIONAL INPUTS)
% 
%      fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed.
%
% OUTPUTS:
%
% geomParameters_flipped: This is a [1XN] vector of the parameters of the
%      flipped geometry
% 
%
% DEPENDENCIES:
%
%   fcn_geometry_plotGeometry
%
% EXAMPLES:
%      
% See the script: script_test_fcn_geometry_flipGeom
% for a full test suite.
%
% This function was written on 2024_05_15 by Sean Brennan
% Questions or comments? sbrennan@psu.edu

% Revision History
% 2024_05_15 - Sean Brennan
% -- wrote the code 

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==3 && isequal(varargin{end},-1))
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

if 0==flag_max_speed
    if flag_check_inputs
        % Are there the right number of inputs?
        narginchk(2,3);

    end
end

% Does user want to specify fig_num?
fig_num = []; % Default is to have no figure
flag_do_plots = 0;
if (0==flag_max_speed) && (3<= nargin)
    temp = varargin{end};
    if ~isempty(temp)
        fig_num = temp;
        flag_do_plots = 1;
    end
end

%% Main Code starts from here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _       
%  |  \/  |     (_)      
%  | \  / | __ _ _ _ __  
%  | |\/| |/ _` | | '_ \ 
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch lower(geomType)
    case 'circle'
        geomParameters_flipped = geomParameters; 
    case 'arc'
        arc_parameters_flipped      = geomParameters;
        arc_parameters_flipped(1,4) = geomParameters(1,5);
        arc_parameters_flipped(1,5) = geomParameters(1,4);
        if 1==geomParameters(1,7)
            arc_parameters_flipped(1,7) = 0;
        else
            arc_parameters_flipped(1,7) = 1;
        end
        geomParameters_flipped = arc_parameters_flipped;

    case 'line'
        line_parameters_flipped         = geomParameters;
        line_parameters_flipped(1,1:2)  = -geomParameters(1,1:2);

        geomParameters_flipped = line_parameters_flipped;
    case {'segment'}
        segment_parameters_flipped         = geomParameters;
        segment_parameters_flipped(1,1:2)  = -geomParameters(1,1:2);
        segment_parameters_flipped(1,5)    = -geomParameters(1,6);
        segment_parameters_flipped(1,6)    = -geomParameters(1,5);

        geomParameters_flipped = segment_parameters_flipped;
    case 'spiral'
        URHERE
        spiral_join_parameters(1,1) = spiralLength;
        spiral_join_parameters(1,2) = -rotation_angle;
        spiral_join_parameters(1,3) = circle1_start_xy(1,1);
        spiral_join_parameters(1,4) = circle1_start_xy(1,2);
        spiral_join_parameters(1,5) = K0;
        spiral_join_parameters(1,6) = Kf;
    otherwise
        warning('on','backtrace');
        warning('An error will be thrown at this point due to missing code.');
        error('fcn_geometry_intersectionGeom is not yet ready for curves from fit type: %s',geomType);

end
%% Plot the results (for debugging)?
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

    temp_h = figure(fig_num);
    flag_rescale_axis = 0;
    if isempty(get(temp_h,'Children'))
        flag_rescale_axis = 1;
    end

    % Plot the inputs
    subplot(1,2,1);
    hold on;
    grid on;
    axis equal;

    title('Inputs')
    xlabel('X [meters]');
    ylabel('Y [meters]')

    fcn_geometry_plotGeometry(lower(geomType),geomParameters,[],sprintf(' ''LineWidth'',4 '));

    % Make axis slightly larger?
    if flag_rescale_axis
        temp = axis;
        %     temp = [min(points(:,1)) max(points(:,1)) min(points(:,2)) max(points(:,2))];
        axis_range_x = temp(2)-temp(1);
        axis_range_y = temp(4)-temp(3);
        percent_larger = 0.3;
        axis([temp(1)-percent_larger*axis_range_x, temp(2)+percent_larger*axis_range_x,  temp(3)-percent_larger*axis_range_y, temp(4)+percent_larger*axis_range_y]);
    end
    temp_axis = axis;

    % Plot the outputs
    subplot(1,2,2);
    hold on;
    grid on;
    axis equal;

    title('Outputs')
    xlabel('X [meters]');
    ylabel('Y [meters]')

    fcn_geometry_plotGeometry(lower(geomType),geomParameters_flipped,[],sprintf(' ''LineWidth'',4 '));

    axis(temp_axis);


    sgtitle(sprintf('Geometry Flipping for: %s',geomType));
    
end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file);
end


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


