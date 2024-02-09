function [breakPointsCell] = fcn_geometry_findBreakpoints(domains)

% test script: script_test_fcn_geometry_findBreakpoints

% This function is still under testing. 
% The instructions will be updated after testing the function

N_houghDomains = size(domains,2) - 1;

% Empty breakPoints domain structure
breakPoints.firstBreakPoint = [nan nan];
breakPoints.lastBreakPoint = [nan nan];
breakPoints.fitType = 'empty';
breakPoints.fitParameters = nan;

% Create a cell array to save the structure of each fit
breakPointsCell = cell(1,N_houghDomains);

% This loop saves the first and last points of each fit domain as the first
% and last break points in breakPointsCell cell array. This cell array also
% stores the fit type and parameters

for i = 1:N_houghDomains

    breakPointsCell{i} = breakPoints;
    breakPointsCell{i}.firstBreakPoint = domains{i}.points_in_domain(1,:);
    breakPointsCell{i}.lastBreakPoint = domains{i}.points_in_domain(end,:);
    breakPointsCell{i}.fitType = domains{i}.best_fit_type;
    breakPointsCell{i}.fitParameters = domains{i}.best_fit_parameters;

end

end