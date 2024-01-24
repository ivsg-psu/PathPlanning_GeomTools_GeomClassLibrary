function [vectorsCloseTOref, SameDirection_vectorsCloseTOref] = fcn_geometry_findVectorsINvicinity(inputVectors, refVector, tolerance)
% fcn_findVectorsINvicinity
%
% Finds the vectors that are in the vicinity of the reference vector and
% are in the same direction of the reference vector
%
% FORMAT: 
% 
% [vectorsCloseTOref, SameDirection_vectorsCloseTOref] = fcn_findVectorsINvicinity(inputVectors, refVector, tolerance)
%
% INPUTS:
%
%   inputVectors: a list of Nxm vector where N is the number of vectors,
%   and m is the dimension of the vector (typically 2 or 3).
%
%   refVector: this is the vector that is used as a reference to check if
%   an input vector is in the vicinity or not.
%
%   tolerance: this is the tolerance limit to determine if a vector is in
%   the vicinity of the reference vector or not. Here, the euclidean
%   distance of the unit vectors of input vectors and the unit vector of
%   reference vector is computed. If the distance is less than tolerance
%   limit, the input vector is considered to be in the vicinity of the
%   reference vector.
% 
% (OPTIONAL INPUTS)
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

flag_do_plots = 0; % Flag to plot
flag_do_debug = 0; % Flag to plot the results for debugging
flag_check_inputs = 1; % Flag to perform input checking

if flag_do_debug
    st = dbstack; %#ok<*UNRCH>
    fprintf(1,'STARTING function: %s, in file: %s\n',st(1).name,st(1).file);
    debug_fig_num = 34838;
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

if flag_check_inputs
    % Are there the right number of inputs?
    narginchk(3,3);

    % Check the points input vectors to be a 2 or 3 column vector
    fcn_DebugTools_checkInputsToFunctions(inputVectors, '2or3column_of_numbers');

    % Check the points reference vector to be a 2 or 3 column vector
    fcn_DebugTools_checkInputsToFunctions(refVector, '2or3column_of_numbers');

    % Check the tolerance input is a positive single number
    fcn_DebugTools_checkInputsToFunctions(tolerance, 'positive_1column_of_numbers',1);

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
unitVectors_inputVectors = (1./sqrt(sum(inputVectors.^2,2))).*inputVectors;

% Calculate the unit vector of the reference vector (refVector)
unitRefVector = refVector./sqrt(sum(refVector.^2,2));

% Calculate the distance between the unit reference vector (unitRefVector)
% and the unit vectors of input vectors (unitVectors_inputVectors)
dist_btw_unit_refVectANDinputVec = sqrt(sum((unitRefVector - unitVectors_inputVectors).^2,2));

% Find indices of inputVectors that are not in the vicinity of refVector
% i.e. distances (dist_btw_unit_refVectANDinputVec) within the tolerance
% limit
% indices = find((dist_btw_unit_refVectANDinputVec <= tolerance) == 1);
indices = (dist_btw_unit_refVectANDinputVec <= tolerance) == 1;

% The inputVectors that are in the vicinity of refVector
vectorsCloseTOref = inputVectors(indices,:);

% The indices of inputVectors that are in the vicinity of refVector and are
% parallel to the refVector
indices_inputVecsPARALLELrefVec = (vectorsCloseTOref./refVector) == vectorsCloseTOref(:,1)./refVector(:,1);
indices_inputVecsPARALLELrefVec = sum(indices_inputVecsPARALLELrefVec,2) == 3; % indicesSameDirection = find(sum(indicesSameDirection,2) == 3);

% The vectors (vectorsCloseTOref) that are in the close vicinity and are
% parallel to the refVector
parallel_vectorsCloseTOref = vectorsCloseTOref(indices_inputVecsPARALLELrefVec,:);

% Determining if the vectors (vectorsCloseTOref) are in the same direction
% to the refVector
checkDirection = (sign(parallel_vectorsCloseTOref) == sign(refVector));

if (checkDirection)
    SameDirection_vectorsCloseTOref = parallel_vectorsCloseTOref;    
else
    indicesSameDirection = (sum(checkDirection,2) == 3);
    SameDirection_vectorsCloseTOref = parallel_vectorsCloseTOref(indicesSameDirection,:);
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





