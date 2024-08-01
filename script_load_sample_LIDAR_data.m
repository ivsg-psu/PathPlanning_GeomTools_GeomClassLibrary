%% script_load_sample_LIDAR_data.m
% Demonstrates the loading of vehicle pose and LIDAR data with many
% conditional checks to confirm the file is present, the variables were not
% already loaded, the files are up-to-date, etc.

% Revision history:
% 2024_07_17 - S Brennan
% -- wrote the code



%% Code inputs
% Set the "inputs" to the file loading process - need the date and names
% and date of file creation for the Vehicle Pose data file

test_date_string = '2024_06_28'; % The date of testing. This defines the folder where the data should be found within LargeData main folder
vehicle_pose_string = 'VehiclePose_ENU.mat'; % The name of the file containing VehiclePose
LIDAR_file_string   = 'Velodyne_LiDAR_Scan_ENU.mat'; % The name of the file containing the LIDAR data
permanent_file_date = '17-Jul-2024 12:38:35'; % The exact date and time of the files

%% Check to see if data was loaded earlier

% Set the file names
% Create the loadfile name for the Vehicle pose data
mat_loadFilename_vehiclePose = fullfile(cd,'LargeData',test_date_string, vehicle_pose_string);
mat_loadFilename_LIDARdata   = fullfile(cd,'LargeData',test_date_string, LIDAR_file_string);


% Set the default flag
flag_load_all_data = 0; % FORCE LOAD? Set this manually to 1 to FORCE load

% If the variable does not yet exist, load it
if ~exist('VehiclePose','var')
    flag_load_all_data = 1;
end
if ~exist('LiDAR_Scan_ENU_Entire_Loop','var')
    flag_load_all_data = 1;
end

% Does the data match in terms of dates?
if exist('permanent_file_date','var') || ~isempty(permanent_file_date)

    % Do the names match? if not, don't use the permanet ones - need to
    % reload!

    % Check the file's date of creation - if it doesn't match, the file has
    % been edited and needs to be reloaded

    % % BELOW is an example, commented out, of how to check date of this
    % % function:
    % st = dbstack;
    % this_function = st(1).file;
    % file_info = dir(which(this_function));

    file_info = dir(which(mat_loadFilename_vehiclePose));
    file_date = file_info.date;

    if ~strcmp(file_date,permanent_file_date)
        fprintf(1,'\n\nComparing the data files creation date to the date given by user, they did not match.\n');
        fprintf(1,'\tCreation date: %s\n',file_date);
        fprintf(1,'\tGiven by user: %s\n',permanent_file_date)
        flag_load_all_data = 1;
    end

else
    flag_load_all_data = 1;
end

%% Load the data?
% If something indiactes that the data is wrong, the flag_load_all_data
% will have been set to 1 in the above steps

if 1==flag_load_all_data

    % Does the file exist?
    if exist(mat_loadFilename_vehiclePose,'file')
        load(mat_loadFilename_vehiclePose,'VehiclePose');
    else
        % File does not exist - need to load it
        error('Unable to find file: %s',mat_loadFilename_vehiclePose);
    end


    % Does the file exist?
    if exist(mat_loadFilename_LIDARdata,'file')
        load(mat_loadFilename_LIDARdata,'LiDAR_Scan_ENU_Entire_Loop');
    else
        % File does not exist - need to load it
        error('Unable to find file: %s',mat_loadFilename_LIDARdata);
    end
else
    fprintf(1,'Vehice pose and LIDAR data appears unchanged. Using existing variables within workspace.\n')
end
