function [orientation, isconvex, Nencirclements] = ...
    fcn_geometry_findPolytopeOrientations(vertex_pts,...
    varargin)
% fcn_geometry_findPolytopeOrientations
% Finds if a 2D polytope is clockwise or counterClockwise and if the
% polytope is convex or non-convex. The number of encirclements is also
% reported, with positive values being positive encirclements and negative
% values being negative encirclements. 
% 
% The method used is to find the direction of the change of orientation of
% each edge to the next, and the angle of change. The direction is
% determined via the cross products of each edge of a polytope with the
% next edge. The angle is determined via the dot product of each edge with
% the next, then taking the arc cosine. The sum of angles is used to
% determine the number of encirclements in the progression around the edges
% of the polytope. The number of encirclements is determined by rounding
% the sum of angle changes to the nearest multiple of 2*pi.
%
% FORMAT:
%
% [orientation, isconvex, Nencirclements] = ...
%     fcn_geometry_findPolytopeOrientations(vertex_pts,...
%     (figNum))
%
% INPUTS:
%
%      vertex_pts: an [N x 2] vector of [X Y] data defining the vertices of
%      the polytope, with N>=3.
%
%      (OPTIONAL INPUTS)
% 
%      figNum: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed.
%
% OUTPUTS:
%
%      orientation: returns 1 if the polytope encloses net
%      counter-clockwise, -1 if the polytope encloses net clockwise
% 
%      isconvex: returns true if the polytope is strictly convex, false otherwise
% 
%      Nencirclements: the number of encirclements in the closure direction 
%
% DEPENDENCIES:
%
%      fcn_DebugTools_checkInputsToFunctions      
%
% EXAMPLES:
%
% See the script: script_test_fcn_geometry_findPolytopeOrientations
% for a full test suite.
%
% This function was written on 2025_10_18 by S. Brennan
% Questions or comments? sbrennan@psu.edu 

% Revision History:
% 2021-10-18
% - First write of the code


%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the figNum variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
MAX_NARGIN = 2; % The largest Number of argument inputs to the function
flag_max_speed = 0;
if (nargin==MAX_NARGIN && isequal(varargin{end},-1))
    flag_do_debug = 0; %   %   Flag to plot the results for debugging
    flag_check_inputs = 0; % Flag to perform input checking
    flag_max_speed = 1;
else
    % Check to see if we are externally setting debug mode to be "on"
    flag_do_debug = 0; %   %   Flag to plot the results for debugging
    flag_check_inputs = 1; % Flag to perform input checking
    MATLABFLAG_GEOMETRY_FLAG_CHECK_INPUTS = getenv("MATLABFLAG_GEOMETRY_FLAG_CHECK_INPUTS");
    MATLABFLAG_GEOMETRY_FLAG_DO_DEBUG = getenv("MATLABFLAG_GEOMETRY_FLAG_DO_DEBUG");
    if ~isempty(MATLABFLAG_GEOMETRY_FLAG_CHECK_INPUTS) && ~isempty(MATLABFLAG_GEOMETRY_FLAG_DO_DEBUG)
        flag_do_debug = str2double(MATLABFLAG_GEOMETRY_FLAG_DO_DEBUG);
        flag_check_inputs  = str2double(MATLABFLAG_GEOMETRY_FLAG_CHECK_INPUTS);
    end
end

% flag_do_debug = 1;

if flag_do_debug
    st = dbstack; %#ok<*UNRCH>
    fprintf(1,'STARTING function: %s, in file: %s\n',st(1).name,st(1).file);
    debug_figNum = 999978;
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
    if 1 == flag_check_inputs

        % Are there the right number of inputs?
        narginchk(1,MAX_NARGIN);

        % Check the vertex_pts input to be at 2 columns, 3+ rows
        fcn_DebugTools_checkInputsToFunctions(...
            vertex_pts, '2column_of_numbers',[3 4]);

    end
end
% Does user want to show the plots?
flag_do_plots = 0; % Default is to show plots
figNum = []; % Empty by default
if (0==flag_max_speed) && (MAX_NARGIN == nargin)
    temp = varargin{end};
    if ~isempty(temp) % Did the user NOT give an empty figure number?
        figNum = temp;
        flag_do_plots = 1;
    end
end


%% Start of main code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _
%  |  \/  |     (_)
%  | \  / | __ _ _ _ __
%  | |\/| |/ _` | | '_ \
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make sure the polytope is closed, e.g. 1st point is equal to last point
if ~isequal(vertex_pts(1,:),vertex_pts(end,:))
    fullCircleEdgePoints = [vertex_pts; vertex_pts(1,:)];
else
    fullCircleEdgePoints = vertex_pts;
end

% For debugging
if flag_do_debug
    figure(debug_figNum); 
    clf;
    hold on;
    axis equal;
    
    plot(fullCircleEdgePoints(:,1),fullCircleEdgePoints(:,2),'r.-');
end

% Take row differences
edgeVectors = diff(fullCircleEdgePoints,1,1);
fullCircleEdgeVectors = [edgeVectors; edgeVectors(1,:)];
magFullCircleEdgeVectors = sum(fullCircleEdgeVectors.^2,2).^0.5;
unitFullCircleEdgeVectors = fullCircleEdgeVectors./magFullCircleEdgeVectors;

% Take dot products
dotProducts = sum(unitFullCircleEdgeVectors(1:end-1,:).*unitFullCircleEdgeVectors(2:end,:),2);
Ncross = length(unitFullCircleEdgeVectors(1:end-1,:));
crossProducts = cross(...
    [unitFullCircleEdgeVectors(1:end-1,:) zeros(Ncross,1)],...
    [unitFullCircleEdgeVectors(2:end,:) zeros(Ncross,1)],2);
crossSigns = sign(crossProducts(:,3));

% Check to see if the polytope is convex
if ~all(crossSigns==crossSigns(1,1))
    isconvex = false;
else
    isconvex = true;
end

angles = real(acos(dotProducts))*180/pi.*crossSigns;

angleSum = sum(angles);

Nencirclements = round(angleSum)/360;

if Nencirclements>=0
    orientation = 1;
else
    orientation = -1;
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
    % check whether the figure already has data
    temp_h = figure(figNum);
    flag_rescale_axis = 0;
    if isempty(get(temp_h,'Children'))
        flag_rescale_axis = 1; 
    end

    hold on;

    magnitudeScale = max(vertex_pts,[],'all') - min(vertex_pts,[],'all');

    % Plot the polytopes
    h_patch = patch(vertex_pts(:,1),vertex_pts(:,2),[0 0 1]);
    set(h_patch,'FaceAlpha',0.7);
    plot(vertex_pts(:,1),vertex_pts(:,2),'b.','MarkerSize',20);

    % Label the points
    addNudge = 0.01*magnitudeScale;
    text(vertex_pts(:,1)+addNudge,vertex_pts(:,2)+addNudge,string((1:length(vertex_pts(:,1)))'));

    titlestring = [];
    if 1==orientation
        titlestring = cat(2,titlestring,'Positive orientation, ');
    else
        titlestring = cat(2,titlestring,'Negative orientation, ');
    end
    if isconvex
        titlestring = cat(2,titlestring,'Convex, ');
    else
        titlestring = cat(2,titlestring,'Non-convex, ');
    end
    titlestring = cat(2,titlestring, sprintf('%.0f',abs(round(Nencirclements))),' Encirclement');
    if Nencirclements~=1
        titlestring = cat(2,titlestring,'s');
    end

    maxY = max(vertex_pts(:,2)) + addNudge*5;
    meanX = mean(vertex_pts(:,1));
    text(meanX, maxY, titlestring,'HorizontalAlignment','center');

    if 1==flag_rescale_axis
        temp = axis;
        axis([temp(1,1)-1 temp(1,2)+1 temp(1,3)-1 temp(1,4)+1]);
    end

end


if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file);
end

end % Ends main function

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

