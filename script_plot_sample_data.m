

%% Check to see if data was loaded earlier
mat_filename = fullfile(cd,'Data','AllMarkerClusterData.mat'); %%%% not loading centerline data
% flag_load_all_data = 0; % Default value

flag_load_all_data = 0; % FORCE LOAD?

if isempty(permanent_MarkerClusterNames) || isempty(permanent_MarkerClusterData) || isempty(permanent_file_date)

    % Does the file exist?
    if exist(mat_filename,'file')
        load(mat_filename,'permanent_MarkerClusterData','permanent_MarkerClusterNames','permanent_file_date');
    else
        % File does not exist - need to load it
        flag_load_all_data = 1;
    end
end

% Does the data match?
if ~isempty(permanent_MarkerClusterNames)

    % Do the names match? if not, don't use the permanet ones - need to
    % reload!
    if ~isequal(MarkerClusterNames,permanent_MarkerClusterNames)
        flag_load_all_data = 1;
    end

    % Check the file's date of creation - if it doesn't match, the file has
    % been edited and needs to be reloaded
    st = dbstack;
    this_function = st(1).file;
    file_info = dir(which(this_function));
    file_date = file_info.date;

    if ~exist('permanent_file_date','var') || ~strcmp(file_date,permanent_file_date)
        flag_load_all_data = 1;
    end

else
    flag_load_all_data = 1;
end
