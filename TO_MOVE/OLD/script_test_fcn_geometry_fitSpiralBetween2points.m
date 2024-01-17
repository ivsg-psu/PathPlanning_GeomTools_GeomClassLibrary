% script_test_fcn_geometry_fitSpiralBetween2points.m
% tests fcn_geometry_fitSpiralBetween2points.m

% Revision history
% 2023_10_19 - Aneesh Batchu
% -- wrote the code originally

%% Set up the workspace

clc
close all

%% Examples for basic path operations and function testing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  ______                           _
% |  ____|                         | |
% | |__  __  ____ _ _ __ ___  _ __ | | ___  ___
% |  __| \ \/ / _` | '_ ` _ \| '_ \| |/ _ \/ __|
% | |____ >  < (_| | | | | | | |_) | |  __/\__ \
% |______/_/\_\__,_|_| |_| |_| .__/|_|\___||___/
%                            | |
%                            |_|
% See: https://patorjk.com/software/taag/#p=display&f=Big&t=Examples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Fit Spiral b/w 2 points - Case 1

% The end point of a fitted straight line 
endPoint_straightLine = [1 2];

% The start point of a fitted arc 
startPoint_arc = [3,4]; 

% The fitted points of an Archimedean Spiral
fittedSpiral = fcn_geometry_fitSpiralBetween2points(endPoint_straightLine, startPoint_arc, 1);
disp(fittedSpiral);

%% Fit Spiral b/w 2 points - Case 2

% The end point of a fitted straight line 
endPoint_straightLine = [1 3];

% The start point of a fitted arc 
startPoint_arc = [3,3]; 

% The fitted points of an Archimedean Spiral
fittedSpiral = fcn_geometry_fitSpiralBetween2points(endPoint_straightLine, startPoint_arc, 12);
disp(fittedSpiral);

%% Fit Spiral b/w 2 points - Case 3

% The end point of a fitted straight line 
endPoint_straightLine = [1 3];

% The start point of a fitted arc 
startPoint_arc = [3,1]; 

% The fitted points of an Archimedean Spiral
fittedSpiral = fcn_geometry_fitSpiralBetween2points(endPoint_straightLine, startPoint_arc, 123);
disp(fittedSpiral);

%% Fit Spiral b/w 2 points - Case 4

% The end point of a fitted straight line 
endPoint_straightLine = [6 1];

% The start point of a fitted arc 
startPoint_arc = [99,56]; 

% The fitted points of an Archimedean Spiral
fittedSpiral = fcn_geometry_fitSpiralBetween2points(endPoint_straightLine, startPoint_arc, 1234);
disp(fittedSpiral);

