function [verifyIntersectionBool, LHs, RHs] = fcn_geometry_verifyIntersectionPoint(endPointsCell, intersectionPoints)

% Need to update the instructions

% 2024_02_29 - Aneesh Batchu
% Wrote the code originally

% Fit parameters of Line Segment
%
% breakPointsCell{}.fitParameters(1) = unit_projection_vector_x 
% breakPointsCell{}.fitParameters(2) = unit_projection_vector_y 
% breakPointsCell{}.fitParameters(3) = base_point_x 
% breakPointsCell{}.fitParameters(4) = base_point_y
% breakPointsCell{}.fitParameters(5) = station_distance_min 
% breakPointsCell{}.fitParameters(6) = station_distance_max

% Creating a matrix to store the best fit parameters
fitParametersMatrix = zeros([length(endPointsCell)-1,6]);

for i = 1: length(endPointsCell)-1
    fitParametersMatrix(i,:) = endPointsCell{i}.fitParameters;
    
end

% Finding the slopes of the lines
slopesLineSegments = fitParametersMatrix(:,2)./fitParametersMatrix(:,1);

% Finding the y-intercepts of the lines 
yInterceptsLines = fitParametersMatrix(:,4) - slopesLineSegments.*fitParametersMatrix(:,3);

RHs =  round(slopesLineSegments.*intersectionPoints(:,1) + yInterceptsLines, 4);
LHs = round(intersectionPoints(:,2),4);

verifyIntersectionBool = (LHs ==  RHs);


end