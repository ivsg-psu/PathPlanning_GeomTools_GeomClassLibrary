%% Introduction to and Purpose of the Geometry Class library
% This is a demonstration script to show the primary functionality of the
% GeomClass library.
%
% This is the explanation of the code that can be found by running
%       script_demo_Geometry.m
% 
% This code repo is typically located at:
%
%   https://github.com/ivsg-psu/PathPlanning_GeomTools_GeomClassLibrary
%
% If you have questions or comments, please contact Sean Brennan at
% sbrennan@psu.edu 
% 
% Additional contributers include:
% 2016 - Seth Tau (now with Army)
% 2023 - Aneesh Batchu



%% Revision History:
% 2023_11_21 - sbrennan@psu.edu
% -- started a demo code set
% 2023_12_27 - sbrennan@psu.edu
% -- switched to using environment settings for checking input parameters

%% Prep the workspace
close all
clc

%% Dependencies and Setup of the Code
% The code requires several other libraries to work, namely the following
%
% * DebugTools - the repo can be found at: https://github.com/ivsg-psu/Errata_Tutorials_DebugTools


% List what libraries we need, and where to find the codes for each
clear library_name library_folders library_url

ith_library = 1;
library_name{ith_library}    = 'DebugTools_v2023_04_22';
library_folders{ith_library} = {'Functions','Data'};
library_url{ith_library}     = 'https://github.com/ivsg-psu/Errata_Tutorials_DebugTools/archive/refs/tags/DebugTools_v2023_04_22.zip';

% ith_library = ith_library+1;
% library_name{ith_library}    = 'PathClass_v2023_10_18';
% library_folders{ith_library} = {'Functions'};
% library_url{ith_library}     = 'https://github.com/ivsg-psu/PathPlanning_PathTools_PathClassLibrary/archive/refs/tags/PathClass_v2023_10_18.zip';
% 
% ith_library = ith_library+1;
% library_name{ith_library}    = 'GPSClass_v2023_06_29';
% library_folders{ith_library} = {'Functions'};
% library_url{ith_library}     = 'https://github.com/ivsg-psu/FieldDataCollection_GPSRelatedCodes_GPSClass/archive/refs/tags/GPSClass_v2023_06_29.zip';
% 
% ith_library = ith_library+1;
% library_name{ith_library}    = 'LineFitting_v2023_07_24';
% library_folders{ith_library} = {'Functions'};
% library_url{ith_library}     = 'https://github.com/ivsg-psu/FeatureExtraction_Association_LineFitting/archive/refs/tags/LineFitting_v2023_07_24.zip';
% 
% ith_library = ith_library+1;
% library_name{ith_library}    = 'FindCircleRadius_v2023_08_02';
% library_folders{ith_library} = {'Functions'};                                
% library_url{ith_library}     = 'https://github.com/ivsg-psu/PathPlanning_GeomTools_FindCircleRadius/archive/refs/tags/FindCircleRadius_v2023_08_02.zip';
% 
% ith_library = ith_library+1;
% library_name{ith_library}    = 'BreakDataIntoLaps_v2023_08_25';
% library_folders{ith_library} = {'Functions'};                                
% library_url{ith_library}     = 'https://github.com/ivsg-psu/FeatureExtraction_DataClean_BreakDataIntoLaps/archive/refs/tags/BreakDataIntoLaps_v2023_08_25.zip';
% 
% ith_library = ith_library+1;
% library_name{ith_library}    = 'ParseXODR_v2023_10_23';
% library_folders{ith_library} = {'Functions'};                                
% library_url{ith_library}     = 'https://github.com/ivsg-psu/PathPlanning_MapTools_ParseXODR/archive/refs/tags/ParseXODR_v2023_10_23.zip';




%% Clear paths and folders, if needed
if 1==0
    clear flag_GeomTools_Folders_Initialized;
    fcn_INTERNAL_clearUtilitiesFromPathAndFolders;

    % Clean up data files
    traces_mat_filename = fullfile(cd,'Data','AllTracesData.mat'); %%%% not loading centerline data
    if exist(traces_mat_filename,'file')
        delete(traces_mat_filename);
    end
    marker_clusters_mat_filename = fullfile(cd,'Data','AllMarkerClusterData.mat'); %%%% not loading centerline data
    if exist(marker_clusters_mat_filename,'file')

        delete(marker_clusters_mat_filename);
    end

end


%% Do we need to set up the work space?
if ~exist('flag_GeomTools_Folders_Initialized','var')
    this_project_folders = {'Functions','Data'};
    fcn_INTERNAL_initializeUtilities(library_name,library_folders,library_url,this_project_folders);
    flag_GeomTools_Folders_Initialized = 1;
end




%% Clear paths and folders, if needed
if 1==0
    clear flag_GeomClass_Folders_Initialized;
    fcn_INTERNAL_clearUtilitiesFromPathAndFolders;
end


%% Do we need to set up the work space?
if ~exist('flag_GeomClass_Folders_Initialized','var')
    this_project_folders = {'Functions','Data'};
    fcn_INTERNAL_initializeUtilities(library_name,library_folders,library_url,this_project_folders);
    flag_GeomClass_Folders_Initialized = 1;
end

%% Set environment flags for input checking
% These are values to set if we want to check inputs or do debugging
% setenv('MATLABFLAG_GEOMETRY_FLAG_CHECK_INPUTS','1');
% setenv('MATLABFLAG_GEOMETRY_FLAG_DO_DEBUG','1');
setenv('MATLABFLAG_GEOMETRY_FLAG_CHECK_INPUTS','1');
setenv('MATLABFLAG_GEOMETRY_FLAG_DO_DEBUG','0');


%% Circle-related calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% See http://patorjk.com/software/taag/#p=display&f=Big&t=Filling%20Test%20Data
 %   _____ _          _                      _       _           _    _____      _            _       _   _                 
 %  / ____(_)        | |                    | |     | |         | |  / ____|    | |          | |     | | (_)                
 % | |     _ _ __ ___| | ___ ______ _ __ ___| | __ _| |_ ___  __| | | |     __ _| | ___ _   _| | __ _| |_ _  ___  _ __  ___ 
 % | |    | | '__/ __| |/ _ \______| '__/ _ \ |/ _` | __/ _ \/ _` | | |    / _` | |/ __| | | | |/ _` | __| |/ _ \| '_ \/ __|
 % | |____| | | | (__| |  __/      | | |  __/ | (_| | ||  __/ (_| | | |___| (_| | | (__| |_| | | (_| | |_| | (_) | | | \__ \
 %  \_____|_|_|  \___|_|\___|      |_|  \___|_|\__,_|\__\___|\__,_|  \_____\__,_|_|\___|\__,_|_|\__,_|\__|_|\___/|_| |_|___/
 % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Calculating a circle from 3 points
fig_num = 102;
points = [0 0; 0.5 4; 1 -1; 4 -3; 6 2; 7 -2; 9 3; 11 3; 15 -0.5];
hold on
figure(1); clf;
for i=1:length(points(:,1))-2
    [centers, radii] = fcn_geometry_circleCenterFrom3Points(points(i:i+2,:),fig_num);
    plot(points(:,1),points(:,2),'r-');
end

%% Plotting a circle
fig_num = 107;
fcn_geometry_plotCircle(centers,radii,'b-',fig_num) 

%% Calculating an arc from 2 points and a direction
fig_num = 100;

centers = [0 0; 4 4; 8 10; -6 10];
radii = [1; 2; 4; 3];
start_angles = [90; 0; -90; 45]*pi/180;
start_points_on_circle = [radii.*cos(start_angles) radii.*sin(start_angles)]+centers;
end_angles = [45; 135; 180; 0]*pi/180;
end_points_on_circle = [radii.*cos(end_angles) radii.*sin(end_angles)]+centers;
cross_products = [-1; 1; -1; 1];

true_angle = start_angles - end_angles;

[angles] = fcn_geometry_findAngleUsing2PointsOnCircle(...
    centers,...
    radii,...
    start_points_on_circle,...
    end_points_on_circle,...
    cross_products,...
    fig_num);

assert(isequal(round(angles,4),round([-pi/4; 3*pi/4; -pi/2; 7*pi/4],4)));

%% Calculating an arc direction from 3 points
fig_num = 102;
points1 = [0 0];
points2 = [-1 4];
points3 = [0 5];
is_counterClockwise = fcn_geometry_arcDirectionFrom3Points(points1, points2, points3, fig_num);

assert(isequal(is_counterClockwise,-1))

%% Calculating an arc from 3 points
fig_num = 103;
Radius = 2;
points = Radius*[[cos(0)    sin(0)]; [cos(pi/4) sin(pi/4)]; [cos(pi/2) sin(pi/2)]];

points1 = points(1,:);
points2 = points(3,:);
points3 = points(2,:);

[arc_angle_in_radians_1_to_2, arc_angle_in_radians_1_to_3, circle_centers, radii, start_angles_in_radians] = ...
    fcn_geometry_arcAngleFrom3Points(points1, points2, points3, fig_num);
assert(isequal(arc_angle_in_radians_1_to_2,-3*pi/2));
assert(isequal(arc_angle_in_radians_1_to_3,-7*pi/4));
assert(isequal(circle_centers,[0 0]));
assert(isequal(radii,2));
assert(isequal(start_angles_in_radians,0));



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% See http://patorjk.com/software/taag/#p=display&f=Big&t=Filling%20Test%20Data
%
%  ______ _ _ _ _               _______        _     _____        _
% |  ____(_) | (_)             |__   __|      | |   |  __ \      | |
% | |__   _| | |_ _ __   __ _     | | ___  ___| |_  | |  | | __ _| |_ __ _
% |  __| | | | | | '_ \ / _` |    | |/ _ \/ __| __| | |  | |/ _` | __/ _` |
% | |    | | | | | | | | (_| |    | |  __/\__ \ |_  | |__| | (_| | || (_| |
% |_|    |_|_|_|_|_| |_|\__, |    |_|\___||___/\__| |_____/ \__,_|\__\__,_|
%                        __/ |
%                       |___/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fig_num = 1;
figure(fig_num);
clf;

seed_points = [2 3; 4 5; 7 0; 9 5];
M = 10;
sigma = 0.02;

test_points = fcn_geometry_fillLineTestPoints(seed_points, M, sigma, fig_num);


% Corrupt the results with outliers
fig_num = 2;
probability_of_corruption = 0.2;
magnitude_of_corruption = 4; % 4 times the y-range

corrupted_test_points = fcn_geometry_corruptPointsWithOutliers(test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (fig_num));

% Shuffle points?
shuffled_corrupted_test_points = fcn_geometry_shufflePointOrdering(corrupted_test_points);

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% See http://patorjk.com/software/taag/#p=display&f=Big&t=Perfect%20Fits
%
%  _____           __          _     ______ _ _
% |  __ \         / _|        | |   |  ____(_) |
% | |__) |__ _ __| |_ ___  ___| |_  | |__   _| |_ ___
% |  ___/ _ \ '__|  _/ _ \/ __| __| |  __| | | __/ __|
% | |  |  __/ |  | ||  __/ (__| |_  | |    | | |_\__ \
% |_|   \___|_|  |_| \___|\___|\__| |_|    |_|\__|___/
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[slope,intercept] = fcn_geometry_fitSlopeInterceptNPoints(seed_points(1:2,:));


%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% See http://patorjk.com/software/taag/#p=display&f=Big&t=Hough%20Fits
 %  _    _                   _       ______ _ _       
 % | |  | |                 | |     |  ____(_) |      
 % | |__| | ___  _   _  __ _| |__   | |__   _| |_ ___ 
 % |  __  |/ _ \| | | |/ _` | '_ \  |  __| | | __/ __|
 % | |  | | (_) | |_| | (_| | | | | | |    | | |_\__ \
 % |_|  |_|\___/ \__,_|\__, |_| |_| |_|    |_|\__|___/
 %                      __/ |                         
 %                     |___/                          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Demo Hough line fitting
fig_num = 1111;
transverse_tolerance = 0.2;
station_tolerance = [];

[fitted_parameters, best_fit_source_indicies, best_agreement_indicies] = ...
    fcn_geometry_fitHoughLine(shuffled_corrupted_test_points, transverse_tolerance, station_tolerance,  fig_num);

%% Demo Hough line segment fitting
fig_num = 11111;
transverse_tolerance = 0.2;
station_tolerance = 0.4;

[fitted_parameters, best_fit_source_indicies, best_agreement_indicies] = ...
    fcn_geometry_fitHoughLine(shuffled_corrupted_test_points, transverse_tolerance, station_tolerance,  fig_num);

%% Demo linear pseudo-regression from Hough Line votes
% Show how to do line fitting from Hough votes
fig_num = 11111; % Reuse the previous figure number, as it puts the results atop that plot

% Extract out points from the indicies
source_points = shuffled_corrupted_test_points(best_fit_source_indicies,:);
associated_points_in_domain = shuffled_corrupted_test_points(best_agreement_indicies,:); 

% Perform the regression fit
[best_fit_parameters, best_fit_domain_box] = ...
    fcn_geometry_calcLinearRegressionFromHoughFit(source_points,associated_points_in_domain, fig_num);

% The following 3 lines show how to convert the domain box into a
% polyshape, and how to query points using isinterior to identify which
% points are inside the regression domain box
domainPolyShape = polyshape(best_fit_domain_box(:,1),best_fit_domain_box(:,2),'Simplify',false,'KeepCollinearPoints',true);
IndiciesOfPointsInDomain = isinterior(domainPolyShape,input_points);
best_fit_associated_indicies = find(IndiciesOfPointsInDomain);

%% Demo Hough circle fitting
seed_points = [2 3; 4 5; 6 3];
[true_circleCenter, true_circleRadius] = fcn_geometry_circleCenterFrom3Points(seed_points(1,:),seed_points(2,:),seed_points(3,:),-1);
trueParameters = [true_circleCenter true_circleRadius];

M = 10; % Number of points per meter
sigma = 0.02;

test_points = fcn_geometry_fillArcTestPoints(seed_points, M, sigma); %, fig_num);

% Add outliers?
% Corrupt the results
probability_of_corruption = 0.3;
magnitude_of_corruption = 1;

test_points = fcn_geometry_corruptPointsWithOutliers(test_points,...
    (probability_of_corruption), (magnitude_of_corruption));

fig_num = 22222;
transverse_tolerance = 0.2;
station_tolerance = 0.4;

[fitted_parameters, agreement_indicies] = fcn_geometry_fitHoughCircle(shuffled_corrupted_test_points, transverse_tolerance, station_tolerance,  fig_num);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% See http://patorjk.com/software/taag/#p=display&f=Big&t=Checking%20Fits
%
 %   _____ _               _    _               ______ _ _       
 %  / ____| |             | |  (_)             |  ____(_) |      
 % | |    | |__   ___  ___| | ___ _ __   __ _  | |__   _| |_ ___ 
 % | |    | '_ \ / _ \/ __| |/ / | '_ \ / _` | |  __| | | __/ __|
 % | |____| | | |  __/ (__|   <| | | | | (_| | | |    | | |_\__ \
 %  \_____|_| |_|\___|\___|_|\_\_|_| |_|\__, | |_|    |_|\__|___/
 %                                       __/ |                   
 %                                      |___/                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



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

%% function fcn_INTERNAL_clearUtilitiesFromPathAndFolders
function fcn_INTERNAL_clearUtilitiesFromPathAndFolders
% Clear out the variables
clear global flag* FLAG*
clear flag*
clear path

% Clear out any path directories under Utilities
path_dirs = regexp(path,'[;]','split');
utilities_dir = fullfile(pwd,filesep,'Utilities');
for ith_dir = 1:length(path_dirs)
    utility_flag = strfind(path_dirs{ith_dir},utilities_dir);
    if ~isempty(utility_flag)
        rmpath(path_dirs{ith_dir});
    end
end

% Delete the Utilities folder, to be extra clean!
if  exist(utilities_dir,'dir')
    [status,message,message_ID] = rmdir(utilities_dir,'s');
    if 0==status
        error('Unable remove directory: %s \nReason message: %s \nand message_ID: %s\n',utilities_dir, message,message_ID);
    end
end

end % Ends fcn_INTERNAL_clearUtilitiesFromPathAndFolders

%% fcn_INTERNAL_initializeUtilities
function  fcn_INTERNAL_initializeUtilities(library_name,library_folders,library_url,this_project_folders)
% Reset all flags for installs to empty
clear global FLAG*

fprintf(1,'Installing utilities necessary for code ...\n');

% Dependencies and Setup of the Code
% This code depends on several other libraries of codes that contain
% commonly used functions. We check to see if these libraries are installed
% into our "Utilities" folder, and if not, we install them and then set a
% flag to not install them again.

% Set up libraries
for ith_library = 1:length(library_name)
    dependency_name = library_name{ith_library};
    dependency_subfolders = library_folders{ith_library};
    dependency_url = library_url{ith_library};
    
    fprintf(1,'\tAdding library: %s ...',dependency_name);
    fcn_INTERNAL_DebugTools_installDependencies(dependency_name, dependency_subfolders, dependency_url);
    clear dependency_name dependency_subfolders dependency_url
    fprintf(1,'Done.\n');
end

% Set dependencies for this project specifically
fcn_DebugTools_addSubdirectoriesToPath(pwd,this_project_folders);

disp('Done setting up libraries, adding each to MATLAB path, and adding current repo folders to path.');
end % Ends fcn_INTERNAL_initializeUtilities


function fcn_INTERNAL_DebugTools_installDependencies(dependency_name, dependency_subfolders, dependency_url, varargin)
%% FCN_DEBUGTOOLS_INSTALLDEPENDENCIES - MATLAB package installer from URL
%
% FCN_DEBUGTOOLS_INSTALLDEPENDENCIES installs code packages that are
% specified by a URL pointing to a zip file into a default local subfolder,
% "Utilities", under the root folder. It also adds either the package
% subfoder or any specified sub-subfolders to the MATLAB path.
%
% If the Utilities folder does not exist, it is created.
%
% If the specified code package folder and all subfolders already exist,
% the package is not installed. Otherwise, the folders are created as
% needed, and the package is installed.
%
% If one does not wish to put these codes in different directories, the
% function can be easily modified with strings specifying the
% desired install location.
%
% For path creation, if the "DebugTools" package is being installed, the
% code installs the package, then shifts temporarily into the package to
% complete the path definitions for MATLAB. If the DebugTools is not
% already installed, an error is thrown as these tools are needed for the
% path creation.
%
% Finally, the code sets a global flag to indicate that the folders are
% initialized so that, in this session, if the code is called again the
% folders will not be installed. This global flag can be overwritten by an
% optional flag input.
%
% FORMAT:
%
%      fcn_DebugTools_installDependencies(...
%           dependency_name, ...
%           dependency_subfolders, ...
%           dependency_url)
%
% INPUTS:
%
%      dependency_name: the name given to the subfolder in the Utilities
%      directory for the package install
%
%      dependency_subfolders: in addition to the package subfoder, a list
%      of any specified sub-subfolders to the MATLAB path. Leave blank to
%      add only the package subfolder to the path. See the example below.
%
%      dependency_url: the URL pointing to the code package.
%
%      (OPTIONAL INPUTS)
%      flag_force_creation: if any value other than zero, forces the
%      install to occur even if the global flag is set.
%
% OUTPUTS:
%
%      (none)
%
% DEPENDENCIES:
%
%      This code will automatically get dependent files from the internet,
%      but of course this requires an internet connection. If the
%      DebugTools are being installed, it does not require any other
%      functions. But for other packages, it uses the following from the
%      DebugTools library: fcn_DebugTools_addSubdirectoriesToPath
%
% EXAMPLES:
%
% % Define the name of subfolder to be created in "Utilities" subfolder
% dependency_name = 'DebugTools_v2023_01_18';
%
% % Define sub-subfolders that are in the code package that also need to be
% % added to the MATLAB path after install; the package install subfolder
% % is NOT added to path. OR: Leave empty ({}) to only add
% % the subfolder path without any sub-subfolder path additions.
% dependency_subfolders = {'Functions','Data'};
%
% % Define a universal resource locator (URL) pointing to the zip file to
% % install. For example, here is the zip file location to the Debugtools
% % package on GitHub:
% dependency_url = 'https://github.com/ivsg-psu/Errata_Tutorials_DebugTools/blob/main/Releases/DebugTools_v2023_01_18.zip?raw=true';
%
% % Call the function to do the install
% fcn_DebugTools_installDependencies(dependency_name, dependency_subfolders, dependency_url)
%
% This function was written on 2023_01_23 by S. Brennan
% Questions or comments? sbrennan@psu.edu

% Revision history:
% 2023_01_23:
% -- wrote the code originally
% 2023_04_20:
% -- improved error handling
% -- fixes nested installs automatically

% TO DO
% -- Add input argument checking

flag_do_debug = 0; % Flag to show the results for debugging
flag_do_plots = 0; % % Flag to plot the final results
flag_check_inputs = 1; % Flag to perform input checking

if flag_do_debug
    st = dbstack; %#ok<*UNRCH>
    fprintf(1,'STARTING function: %s, in file: %s\n',st(1).name,st(1).file);
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
    narginchk(3,4);
end

%% Set the global variable - need this for input checking
% Create a variable name for our flag. Stylistically, global variables are
% usually all caps.
flag_varname = upper(cat(2,'flag_',dependency_name,'_Folders_Initialized'));

% Make the variable global
eval(sprintf('global %s',flag_varname));

if nargin==4
    if varargin{1}
        eval(sprintf('clear global %s',flag_varname));
    end
end

%% Main code starts here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _
%  |  \/  |     (_)
%  | \  / | __ _ _ _ __
%  | |\/| |/ _` | | '_ \
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



if ~exist(flag_varname,'var') || isempty(eval(flag_varname))
    % Save the root directory, so we can get back to it after some of the
    % operations below. We use the Print Working Directory command (pwd) to
    % do this. Note: this command is from Unix/Linux world, but is so
    % useful that MATLAB made their own!
    root_directory_name = pwd;
    
    % Does the directory "Utilities" exist?
    utilities_folder_name = fullfile(root_directory_name,'Utilities');
    if ~exist(utilities_folder_name,'dir')
        % If we are in here, the directory does not exist. So create it
        % using mkdir
        [success_flag,error_message,message_ID] = mkdir(root_directory_name,'Utilities');
        
        % Did it work?
        if ~success_flag
            error('Unable to make the Utilities directory. Reason: %s with message ID: %s\n',error_message,message_ID);
        elseif ~isempty(error_message)
            warning('The Utilities directory was created, but with a warning: %s\n and message ID: %s\n(continuing)\n',error_message, message_ID);
        end
        
    end
    
    % Does the directory for the dependency folder exist?
    dependency_folder_name = fullfile(root_directory_name,'Utilities',dependency_name);
    if ~exist(dependency_folder_name,'dir')
        % If we are in here, the directory does not exist. So create it
        % using mkdir
        [success_flag,error_message,message_ID] = mkdir(utilities_folder_name,dependency_name);
        
        % Did it work?
        if ~success_flag
            error('Unable to make the dependency directory: %s. Reason: %s with message ID: %s\n',dependency_name, error_message,message_ID);
        elseif ~isempty(error_message)
            warning('The %s directory was created, but with a warning: %s\n and message ID: %s\n(continuing)\n',dependency_name, error_message, message_ID);
        end
        
    end
    
    % Do the subfolders exist?
    flag_allFoldersThere = 1;
    if isempty(dependency_subfolders{1})
        flag_allFoldersThere = 0;
    else
        for ith_folder = 1:length(dependency_subfolders)
            subfolder_name = dependency_subfolders{ith_folder};
            
            % Create the entire path
            subfunction_folder = fullfile(root_directory_name, 'Utilities', dependency_name,subfolder_name);
            
            % Check if the folder and file exists that is typically created when
            % unzipping.
            if ~exist(subfunction_folder,'dir')
                flag_allFoldersThere = 0;
            end
        end
    end
    
    % Do we need to unzip the files?
    if flag_allFoldersThere==0
        % Files do not exist yet - try unzipping them.
        save_file_name = tempname(root_directory_name);
        zip_file_name = websave(save_file_name,dependency_url);
        % CANT GET THIS TO WORK --> unzip(zip_file_url, debugTools_folder_name);
        
        % Is the file there?
        if ~exist(zip_file_name,'file')
            error(['The zip file: %s for dependency: %s did not download correctly.\n' ...
                'This is usually because permissions are restricted on ' ...
                'the current directory. Check the code install ' ...
                '(see README.md) and try again.\n'],zip_file_name, dependency_name);
        end
        
        % Try unzipping
        unzip(zip_file_name, dependency_folder_name);
        
        % Did this work? If so, directory should not be empty
        directory_contents = dir(dependency_folder_name);
        if isempty(directory_contents)
            error(['The necessary dependency: %s has an error in install ' ...
                'where the zip file downloaded correctly, ' ...
                'but the unzip operation did not put any content ' ...
                'into the correct folder. ' ...
                'This suggests a bad zip file or permissions error ' ...
                'on the local computer.\n'],dependency_name);
        end
        
        % Check if is a nested install (for example, installing a folder
        % "Toolsets" under a folder called "Toolsets"). This can be found
        % if there's a folder whose name contains the dependency_name
        flag_is_nested_install = 0;
        for ith_entry = 1:length(directory_contents)
            if contains(directory_contents(ith_entry).name,dependency_name)
                if directory_contents(ith_entry).isdir
                    flag_is_nested_install = 1;
                    install_directory_from = fullfile(directory_contents(ith_entry).folder,directory_contents(ith_entry).name);
                    install_files_from = fullfile(directory_contents(ith_entry).folder,directory_contents(ith_entry).name,'*'); % BUG FIX - For Macs, must be *, not *.*
                    install_location_to = fullfile(directory_contents(ith_entry).folder);
                end
            end
        end
        
        if flag_is_nested_install
            [status,message,message_ID] = movefile(install_files_from,install_location_to);
            if 0==status
                error(['Unable to move files from directory: %s\n ' ...
                    'To: %s \n' ...
                    'Reason message: %s\n' ...
                    'And message_ID: %s\n'],install_files_from,install_location_to, message,message_ID);
            end
            [status,message,message_ID] = rmdir(install_directory_from);
            if 0==status
                error(['Unable remove directory: %s \n' ...
                    'Reason message: %s \n' ...
                    'And message_ID: %s\n'],install_directory_from,message,message_ID);
            end
        end
        
        % Make sure the subfolders were created
        flag_allFoldersThere = 1;
        if ~isempty(dependency_subfolders{1})
            for ith_folder = 1:length(dependency_subfolders)
                subfolder_name = dependency_subfolders{ith_folder};
                
                % Create the entire path
                subfunction_folder = fullfile(root_directory_name, 'Utilities', dependency_name,subfolder_name);
                
                % Check if the folder and file exists that is typically created when
                % unzipping.
                if ~exist(subfunction_folder,'dir')
                    flag_allFoldersThere = 0;
                end
            end
        end
        % If any are not there, then throw an error
        if flag_allFoldersThere==0
            error(['The necessary dependency: %s has an error in install, ' ...
                'or error performing an unzip operation. The subfolders ' ...
                'requested by the code were not found after the unzip ' ...
                'operation. This suggests a bad zip file, or a permissions ' ...
                'error on the local computer, or that folders are ' ...
                'specified that are not present on the remote code ' ...
                'repository.\n'],dependency_name);
        else
            % Clean up the zip file
            delete(zip_file_name);
        end
        
    end
    
    
    % For path creation, if the "DebugTools" package is being installed, the
    % code installs the package, then shifts temporarily into the package to
    % complete the path definitions for MATLAB. If the DebugTools is not
    % already installed, an error is thrown as these tools are needed for the
    % path creation.
    %
    % In other words: DebugTools is a special case because folders not
    % added yet, and we use DebugTools for adding the other directories
    if strcmp(dependency_name(1:10),'DebugTools')
        debugTools_function_folder = fullfile(root_directory_name, 'Utilities', dependency_name,'Functions');
        
        % Move into the folder, run the function, and move back
        cd(debugTools_function_folder);
        fcn_DebugTools_addSubdirectoriesToPath(dependency_folder_name,dependency_subfolders);
        cd(root_directory_name);
    else
        try
            fcn_DebugTools_addSubdirectoriesToPath(dependency_folder_name,dependency_subfolders);
        catch
            error(['Package installer requires DebugTools package to be ' ...
                'installed first. Please install that before ' ...
                'installing this package']);
        end
    end
    
    
    % Finally, the code sets a global flag to indicate that the folders are
    % initialized.  Check this using a command "exist", which takes a
    % character string (the name inside the '' marks, and a type string -
    % in this case 'var') and checks if a variable ('var') exists in matlab
    % that has the same name as the string. The ~ in front of exist says to
    % do the opposite. So the following command basically means: if the
    % variable named 'flag_CodeX_Folders_Initialized' does NOT exist in the
    % workspace, run the code in the if statement. If we look at the bottom
    % of the if statement, we fill in that variable. That way, the next
    % time the code is run - assuming the if statement ran to the end -
    % this section of code will NOT be run twice.
    
    eval(sprintf('%s = 1;',flag_varname));
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
    
    % Nothing to do!
    
    
    
end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file);
end

end % Ends function fcn_DebugTools_installDependencies

