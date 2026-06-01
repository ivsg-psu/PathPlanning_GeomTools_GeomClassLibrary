% script_test_fcn_geometry_checkIfPointOn3DPatch
% Exercises the function: fcn_geometry_checkIfPointOn3DPatch
% Revision history:
% 2021_06_12
% - wrote the code
% - revised from fcn_geometry_fitPlaneLinearRegression

close all;


%% Demonstration Examples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  _____                                 _             _   _               ______                           _
% |  __ \                               | |           | | (_)             |  ____|                         | |
% | |  | | ___ _ __ ___   ___  _ __  ___| |_ _ __ __ _| |_ _  ___  _ __   | |__  __  ____ _ _ __ ___  _ __ | | ___  ___
% | |  | |/ _ \ '_ ` _ \ / _ \| '_ \/ __| __| '__/ _` | __| |/ _ \| '_ \  |  __| \ \/ / _` | '_ ` _ \| '_ \| |/ _ \/ __|
% | |__| |  __/ | | | | | (_) | | | \__ \ |_| | | (_| | |_| | (_) | | | | | |____ >  < (_| | | | | | | |_) | |  __/\__ \
% |_____/ \___|_| |_| |_|\___/|_| |_|___/\__|_|  \__,_|\__|_|\___/|_| |_| |______/_/\_\__,_|_| |_| |_| .__/|_|\___||___/
%                                                                                                    | |
%                                                                                                    |_|
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Demonstration%20Examples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Demo case 1: basic plane test in X=Z
figNum = 0001;
figure(figNum);
clf;

patchPoints = [0 0 0; 4 0 4; 3 5 3; 2 2 2; 1 6 1];

NtestPoints = 200;
testPointsX = rand(NtestPoints,1)*10 - 2.5*ones(NtestPoints,1);  
testPointsY = rand(NtestPoints,1)*10 - 2.5*ones(NtestPoints,1) ;  
testPointsZ = testPointsX + randn(NtestPoints,1);
testPoints = [testPointsX testPointsY testPointsZ];

tolerance = 1; 

[flag_isInside3Dpatch, flag_isOn3DPatchPlane, flag_projectsInsidePatch] = fcn_geometry_checkIfPointOn3DPatch(patchPoints, testPoints, (tolerance), (figNum));


% Check variable types
assert(islogical(flag_isInside3Dpatch));
assert(islogical(flag_isOn3DPatchPlane));
assert(islogical(flag_projectsInsidePatch));

% Check variable lengths
Npoints = length(testPoints(:,1));
assert(isequal(size(flag_isInside3Dpatch),[Npoints 1]));
assert(isequal(size(flag_isOn3DPatchPlane),[Npoints 1]));
assert(isequal(size(flag_projectsInsidePatch),[Npoints 1]));

% % Check variable values
% assert(isequal(flag_isInside3Dpatch,[1 0 0 1 0]'));
% assert(isequal(flag_isOn3DPatchPlane,[1 1 0 1 0]'));
% assert(isequal(flag_projectsInsidePatch,[1 0 1 1 0]'));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),figNum));


%% Basic testing examples in 3D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  ____            _        _______        _   _                ______                           _                       ____  _____
% |  _ \          (_)      |__   __|      | | (_)              |  ____|                         | |                     |___ \|  __ \
% | |_) | __ _ ___ _  ___     | | ___  ___| |_ _ _ __   __ _   | |__  __  ____ _ _ __ ___  _ __ | | ___  ___    ______    __) | |  | |
% |  _ < / _` / __| |/ __|    | |/ _ \/ __| __| | '_ \ / _` |  |  __| \ \/ / _` | '_ ` _ \| '_ \| |/ _ \/ __|  |______|  |__ <| |  | |
% | |_) | (_| \__ \ | (__     | |  __/\__ \ |_| | | | | (_| |  | |____ >  < (_| | | | | | | |_) | |  __/\__ \            ___) | |__| |
% |____/ \__,_|___/_|\___|    |_|\___||___/\__|_|_| |_|\__, |  |______/_/\_\__,_|_| |_| |_| .__/|_|\___||___/           |____/|_____/
%                                                       __/ |                             | |
%                                                      |___/                              |_|
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Basic%20Testing%20%20Examples%20%20-%203D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Test case 1: basic plane test in XY
figNum = 1001;
figure(figNum);
clf;

patchPoints = [0 0 0; 10 0 0; 10 10 0]*0.5;

testPoints = [...
    2 1 0; % In the plane and in the enclosed space
    -3 0 0; % In the plane, not enclosed
    5 4 3; % Not in the plane, but enclosed
    9 4 0.1; % Not in the plane, but enclosed and close to plane
    12 12 1]*0.5;  % Not in the plane, and not enclosed

tolerance = []; % Use default

[flag_isInside3Dpatch, flag_isOn3DPatchPlane, flag_projectsInsidePatch] = fcn_geometry_checkIfPointOn3DPatch(patchPoints, testPoints, (tolerance), (figNum));


% Check variable types
assert(islogical(flag_isInside3Dpatch));
assert(islogical(flag_isOn3DPatchPlane));
assert(islogical(flag_projectsInsidePatch));

% Check variable lengths
Npoints = length(testPoints(:,1));
assert(isequal(size(flag_isInside3Dpatch),[Npoints 1]));
assert(isequal(size(flag_isOn3DPatchPlane),[Npoints 1]));
assert(isequal(size(flag_projectsInsidePatch),[Npoints 1]));

% Check variable values
assert(isequal(flag_isInside3Dpatch,[1 0 0 0 0]'));
assert(isequal(flag_isOn3DPatchPlane,[1 1 0 0 0]'));
assert(isequal(flag_projectsInsidePatch,[1 0 1 1 0]'));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),figNum));

%% Test case 2: basic plane test in XY with tolerance
figNum = 1002;
figure(figNum);
clf;

patchPoints = [0 0 0; 10 0 0; 10 10 0]*0.5;

testPoints = [...
    2 1 0; % In the plane and in the enclosed space
    -3 0 0; % In the plane, not enclosed
    5 4 3; % Not in the plane, but enclosed
    9 4 1; % Not in the plane, but enclosed and close to plane (within tolerance)
    12 12 2]*0.5;  % Not in the plane, and not enclosed

tolerance = 0.6; 

[flag_isInside3Dpatch, flag_isOn3DPatchPlane, flag_projectsInsidePatch] = fcn_geometry_checkIfPointOn3DPatch(patchPoints, testPoints, (tolerance), (figNum));


% Check variable types
assert(islogical(flag_isInside3Dpatch));
assert(islogical(flag_isOn3DPatchPlane));
assert(islogical(flag_projectsInsidePatch));

% Check variable lengths
Npoints = length(testPoints(:,1));
assert(isequal(size(flag_isInside3Dpatch),[Npoints 1]));
assert(isequal(size(flag_isOn3DPatchPlane),[Npoints 1]));
assert(isequal(size(flag_projectsInsidePatch),[Npoints 1]));

% Check variable values
assert(isequal(flag_isInside3Dpatch,[1 0 0 1 0]'));
assert(isequal(flag_isOn3DPatchPlane,[1 1 0 1 0]'));
assert(isequal(flag_projectsInsidePatch,[1 0 1 1 0]'));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),figNum));

%% Test case 3: basic plane test in XZ
figNum = 1003;
figure(figNum);
clf;

patchPoints = [0 0 0; 10 0 0; 10 10 0]*0.5;

testPoints = [...
    2 1 0; % In the plane and in the enclosed space
    -3 0 0; % In the plane, not enclosed
    5 4 3; % Not in the plane, but enclosed
    9 4 0.1; % Not in the plane, but enclosed and close to plane
    12 12 1]*0.5;  % Not in the plane, and not enclosed

patchPoints = [patchPoints(:,1) patchPoints(:,3) patchPoints(:,2)];
testPoints  = [testPoints(:,1)  testPoints(:,3)  testPoints(:,2)];
 
tolerance = []; % Use default

[flag_isInside3Dpatch, flag_isOn3DPatchPlane, flag_projectsInsidePatch] = fcn_geometry_checkIfPointOn3DPatch(patchPoints, testPoints, (tolerance), (figNum));


% Check variable types
assert(islogical(flag_isInside3Dpatch));
assert(islogical(flag_isOn3DPatchPlane));
assert(islogical(flag_projectsInsidePatch));

% Check variable lengths
Npoints = length(testPoints(:,1));
assert(isequal(size(flag_isInside3Dpatch),[Npoints 1]));
assert(isequal(size(flag_isOn3DPatchPlane),[Npoints 1]));
assert(isequal(size(flag_projectsInsidePatch),[Npoints 1]));

% Check variable values
assert(isequal(flag_isInside3Dpatch,[1 0 0 0 0]'));
assert(isequal(flag_isOn3DPatchPlane,[1 1 0 0 0]'));
assert(isequal(flag_projectsInsidePatch,[1 0 1 1 0]'));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),figNum));

%% Fast Mode Tests
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  ______        _     __  __           _        _______        _
% |  ____|      | |   |  \/  |         | |      |__   __|      | |
% | |__ __ _ ___| |_  | \  / | ___   __| | ___     | | ___  ___| |_ ___
% |  __/ _` / __| __| | |\/| |/ _ \ / _` |/ _ \    | |/ _ \/ __| __/ __|
% | | | (_| \__ \ |_  | |  | | (_) | (_| |  __/    | |  __/\__ \ |_\__ \
% |_|  \__,_|___/\__| |_|  |_|\___/ \__,_|\___|    |_|\___||___/\__|___/
%
%
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Fast%20Mode%20Tests
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Basic example - NO FIGURE
figNum = 9901;
figure(figNum);
close(figNum);

patchPoints = [0 0 0; 4 0 4; 3 5 3; 2 2 2; 1 6 1];

NtestPoints = 200;
testPointsX = rand(NtestPoints,1)*10 - 2.5*ones(NtestPoints,1);  
testPointsY = rand(NtestPoints,1)*10 - 2.5*ones(NtestPoints,1) ;  
testPointsZ = testPointsX + randn(NtestPoints,1);
testPoints = [testPointsX testPointsY testPointsZ];

tolerance = 1; 

[flag_isInside3Dpatch, flag_isOn3DPatchPlane, flag_projectsInsidePatch] = ...
    fcn_geometry_checkIfPointOn3DPatch(patchPoints, testPoints, (tolerance), ([]));


% Check variable types
assert(islogical(flag_isInside3Dpatch));
assert(islogical(flag_isOn3DPatchPlane));
assert(islogical(flag_projectsInsidePatch));

% Check variable lengths
Npoints = length(testPoints(:,1));
assert(isequal(size(flag_isInside3Dpatch),[Npoints 1]));
assert(isequal(size(flag_isOn3DPatchPlane),[Npoints 1]));
assert(isequal(size(flag_projectsInsidePatch),[Npoints 1]));


% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==figNum));

%% Basic example - NO FIGURE, FAST MODE
figNum = 9902;
figure(figNum);
close(figNum);

patchPoints = [0 0 0; 4 0 4; 3 5 3; 2 2 2; 1 6 1];

NtestPoints = 200;
testPointsX = rand(NtestPoints,1)*10 - 2.5*ones(NtestPoints,1);  
testPointsY = rand(NtestPoints,1)*10 - 2.5*ones(NtestPoints,1) ;  
testPointsZ = testPointsX + randn(NtestPoints,1);
testPoints = [testPointsX testPointsY testPointsZ];

tolerance = 1; 

[flag_isInside3Dpatch, flag_isOn3DPatchPlane, flag_projectsInsidePatch] = ...
    fcn_geometry_checkIfPointOn3DPatch(patchPoints, testPoints, (tolerance), (-1));


% Check variable types
assert(islogical(flag_isInside3Dpatch));
assert(islogical(flag_isOn3DPatchPlane));
assert(islogical(flag_projectsInsidePatch));

% Check variable lengths
Npoints = length(testPoints(:,1));
assert(isequal(size(flag_isInside3Dpatch),[Npoints 1]));
assert(isequal(size(flag_isOn3DPatchPlane),[Npoints 1]));
assert(isequal(size(flag_projectsInsidePatch),[Npoints 1]));

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==figNum));

%% Compare speeds of pre-calculation versus post-calculation versus a fast variant
figNum = 9903;
figure(figNum);
close(figNum);
rng(1823);


patchPoints = [0 0 0; 4 0 4; 3 5 3; 2 2 2; 1 6 1];

NtestPoints = 200;
testPointsX = rand(NtestPoints,1)*10 - 2.5*ones(NtestPoints,1);  
testPointsY = rand(NtestPoints,1)*10 - 2.5*ones(NtestPoints,1) ;  
testPointsZ = testPointsX + randn(NtestPoints,1);
testPoints = [testPointsX testPointsY testPointsZ];

tolerance = 1; 


% Perform the calculation in slow mode
figNum = [];
Niterations = 100; minTimeSlow = Inf; 
tic;
for i=1:Niterations
    tstart = tic;
    [flag_isInside3Dpatch, flag_isOn3DPatchPlane, flag_projectsInsidePatch] = ...
        fcn_geometry_checkIfPointOn3DPatch(patchPoints, testPoints, (tolerance), (figNum));

    telapsed = toc(tstart);
    minTimeSlow = min(telapsed,minTimeSlow);
end
slow_method = toc/Niterations;

% Perform the operation in fast mode
figNum = -1;
minTimeFast = Inf; nsum = 10;
tic;
for i=1:Niterations
    tstart = tic;
    [flag_isInside3Dpatch, flag_isOn3DPatchPlane, flag_projectsInsidePatch] = ...
        fcn_geometry_checkIfPointOn3DPatch(patchPoints, testPoints, (tolerance), (figNum));
    telapsed = toc(tstart);
    minTimeFast = min(telapsed,minTimeFast);
end
fast_method = toc/Niterations;

fprintf(1,'\n\nComparison of fast and slow modes of fcn_geometry_checkIfPointOn3DPatch:\n');
fprintf(1,'N repetitions: %.0d\n',Niterations);
fprintf(1,'Slow mode average speed per call (seconds): %.5f\n',slow_method);
fprintf(1,'Slow mode fastest speed over all calls (seconds): %.5f\n',minTimeSlow);
fprintf(1,'Fast mode average speed per call (seconds): %.5f\n',fast_method);
fprintf(1,'Fast mode fastest speed over all calls (seconds): %.5f\n',minTimeFast);
fprintf(1,'Average ratio of fast mode to slow mode (unitless): %.3f\n',slow_method/fast_method);
fprintf(1,'Fastest ratio of fast mode to slow mode (unitless): %.3f\n',minTimeSlow/minTimeFast);



% Plot results as bar chart
figure(373737);
clf;
X = categorical({'Normal mode','Fast mode'});
X = reordercats(X,{'Normal mode','Fast mode'}); % Forces bars to appear in this exact order, not alphabetized
Y = [slow_method fast_method ]*1000/Niterations;
bar(X,Y)
ylabel('Execution time (Milliseconds)')

%% Fail conditions
if 1==0
    %% FAIL 1: points not long enough
    points = [2 3];
    [root_point, unit_vector] = fcn_geometry_checkIfPointOn3DPatch(points,figNum);
    fprintf(1,'\n\nRoot point is: %.2f %.2f, Unit vector is: %.2f %.2f\n',root_point(1,1),root_point(1,2),unit_vector(1,1),unit_vector(1,2));
end
