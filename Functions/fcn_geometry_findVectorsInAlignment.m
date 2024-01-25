function [dist_btw_unit_refVectAndInputVec, vectorsCloseToRef] = fcn_geometry_findVectorsInAlignment(refVector, inputVectors, tolerance, varargin)
% fcn_findVectorsINvicinity
%
% Finds the vectors that are in alignment with the reference vector. The
% unit vectors of the input vectors and reference are calculated. Then, the
% euclidean distance between the unit vectors of input vectors and unit
% vector of reference vector is calcluated. If the distances are within
% (less than or equal to) the tolerance limit, the corresponding input
% vectors are considered to be aligned with the reference vector
%
% FORMAT: 
% 
% [dist_btw_unit_refVectAndInputVec, vectorsCloseTOref] = fcn_findVectorsINvicinity(inputVectors, refVector, tolerance, fig_num)
%
% INPUTS:
%
%   refVector: this is the vector that is used as a reference to check if
%   an input vector is in alignment or not.
%
%   inputVectors: a list of Nxm vector where N is the number of vectors,
%   and m is the dimension of the vector (typically 2 or 3). Can be of any
%   size.
%   
%   tolerance: this is the tolerance limit to determine if a vector is in
%   alignment with the reference vector or not. Here, the euclidean
%   distance of the unit vectors of input vectors and the unit vector of
%   reference vector is computed. If the distance is less than or equal to
%   tolerance limit, the input vector is considered to be in the vicinity
%   of the reference vector.
% 
% (OPTIONAL INPUTS)
% 
%      fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed.
% 
%      NOTE: Figure can only be plotted if the reference vector and input
%      vectors have exactly three columns. 
%
% OUTPUTS: 
%  
%   vectorsCloseTOref: the vectors that are in the vicinity of the
%   reference vector are calculated
%
%   SameDirection_vectorsCloseTOref: the vectors that are in the vicinity
%   and are in the same direction of the reference vector are computed
%   
%
% DEPENDENCIES:
%
%   fcn_DebugTools_checkInputsToFunctions
%
% EXAMPLES:
%      
% See the script: script_test_fcn_findVectorsINvicinity
% for a full test suite.
%
% This function was written on 2024_01_23 by Aneesh Batchu
% Questions or comments? abb6486@psu.edu or sbrennan@psu.edu

% Revision history:
% 2024_01_23 
% -- wrote the code - Aneesh Batchu

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
        narginchk(3,4);
        
        % Test if reference vector and input Vectors are numeric inputs
        if ~isnumeric(refVector) || ~isnumeric(inputVectors)
            error('The refVector and inputVectors must be numeric inputs')
        end

        if length(refVector(1,:)) ~= length(inputVectors(1,:))
            error('The number of columns in inputVectors must be the same as the number of columns in refVectors.')
        end

        % Check the tolerance input is a positive single number
        fcn_DebugTools_checkInputsToFunctions(tolerance, 'positive_1column_of_numbers',1);

    end
end

% Does user want to specify fig_num?
fig_num = []; % Default is to have no figure
flag_do_plots = 0;
if (0==flag_max_speed) && (4<= nargin)
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

% Calculate the unit vectors of the input vectors (inputVectors)
% unitVectors_inputVectors = fcn_geometry_calcUnitVector(input_vectors, fig_num);
% unitVectors_inputVectors = (1./sqrt(sum(inputVectors.^2,2))).*inputVectors;
vector_magnitude = sum(inputVectors.^2,2).^0.5;
unitVectors_inputVectors = inputVectors./vector_magnitude;

% Calculate the unit vector of the reference vector (refVector)
unitRefVector = refVector./sqrt(sum(refVector.^2,2));

% Calculate the distance between the unit reference vector (unitRefVector)
% and the unit vectors of input vectors (unitVectors_inputVectors)
% Make this one of the inputs
dist_btw_unit_refVectAndInputVec = sqrt(sum((unitRefVector - unitVectors_inputVectors).^2,2));

% Find indices of inputVectors that are not in the vicinity of refVector
% i.e. distances (dist_btw_unit_refVectANDinputVec) within the tolerance
% limit
indices = (dist_btw_unit_refVectAndInputVec <= tolerance) == 1;

% The inputVectors that are in the vicinity of refVector
vectorsCloseToRef = inputVectors(indices,:);

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

    figure(fig_num)
    hold on;
    grid on;
    % axis equal
    
    % origin for plotting
    origin = [0 0 0];

    %plot the unit reference vector in GREEN
    quiver3(origin(1,1), origin(1,2), origin(1,3), unitRefVector(1,1), unitRefVector(1,2), unitRefVector(1,3), 'Color', 'g', 'LineWidth', 4);

    % plot the unit vectors of the input vectors in RED
    N_inputVectors = length(unitVectors_inputVectors(:,1));

    for ith_vector = 1:N_inputVectors

        quiver3(origin(1,1), origin(1,2), origin(1,3), unitVectors_inputVectors(ith_vector,1), unitVectors_inputVectors(ith_vector,2), unitVectors_inputVectors(ith_vector,3), 'Color', 'r', 'LineWidth', 1);

    end

    % plot the vectors that are in alignment with reference vector in BLUE
    N_vectorsCloseTOref = length(vectorsCloseToRef(:,1));

    for ith_vector = 1:N_vectorsCloseTOref

        quiver3(origin(1,1), origin(1,2), origin(1,3), vectorsCloseToRef(ith_vector,1), vectorsCloseToRef(ith_vector,2), vectorsCloseToRef(ith_vector,3), 'Color', 'b', 'LineWidth', 2);

    end

    xlabel('X-axis');
    ylabel('Y-axis');
    zlabel('Z-axis');
    
    view(3) % View plot in 3D

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





