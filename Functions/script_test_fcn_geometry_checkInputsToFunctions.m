% script_test_fcn_geometry_checkInputsToFunctions.m
% Tests fcn_geometry_checkInputsToFunctions

% Revision history:
%      2021_04_22:
%      -- first write of the code copying functionality from fcn_FastestTraversal_checkInputsToFunctions

%% column_of_numbers
%             _                               __                       _
%            | |                             / _|                     | |
%    ___ ___ | |_   _ _ __ ___  _ __    ___ | |_ _ __  _   _ _ __ ___ | |__   ___ _ __ ___
%   / __/ _ \| | | | | '_ ` _ \| '_ \  / _ \|  _| '_ \| | | | '_ ` _ \| '_ \ / _ \ '__/ __|
%  | (_| (_) | | |_| | | | | | | | | || (_) | | | | | | |_| | | | | | | |_) |  __/ |  \__ \
%   \___\___/|_|\__,_|_| |_| |_|_| |_| \___/|_| |_| |_|\__,_|_| |_| |_|_.__/ \___|_|  |___/
%                                  ______   ______
%                                 |______| |______|
%
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Inputs

%% Test the column_of_numbers type (success)
column_of_numbers_test = 4;
fcn_geometry_checkInputsToFunctions(column_of_numbers_test, 'column_of_numbers');

column_of_numbers_test = [4; 3; 2];
fcn_geometry_checkInputsToFunctions(column_of_numbers_test, 'column_of_numbers');

column_of_numbers_test = [4; 3; 2];
fcn_geometry_checkInputsToFunctions(column_of_numbers_test, 'column_of_numbers',3);

assert(true); % pass the test defined by this section if no errors were thrown

%% 2column_of_numbers
% 
%   ___           _                               __                       _                   
%  |__ \         | |                             / _|                     | |                  
%     ) |___ ___ | |_   _ _ __ ___  _ __    ___ | |_ _ __  _   _ _ __ ___ | |__   ___ _ __ ___ 
%    / // __/ _ \| | | | | '_ ` _ \| '_ \  / _ \|  _| '_ \| | | | '_ ` _ \| '_ \ / _ \ '__/ __|
%   / /| (_| (_) | | |_| | | | | | | | | || (_) | | | | | | |_| | | | | | | |_) |  __/ |  \__ \
%  |____\___\___/|_|\__,_|_| |_| |_|_| |_| \___/|_| |_| |_|\__,_|_| |_| |_|_.__/ \___|_|  |___/
%                                      ______   ______                                         
%                                     |______| |______|                                        
%
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Inputs

%% Test the column_of_numbers type (success)
Twocolumn_of_numbers_test = [4 2];
fcn_geometry_checkInputsToFunctions(Twocolumn_of_numbers_test, '2column_of_numbers');

Twocolumn_of_numbers_test = [4 1; 3 0; 2 5];
fcn_geometry_checkInputsToFunctions(Twocolumn_of_numbers_test, '2column_of_numbers');

Twocolumn_of_numbers_test = [4 1; 3 9; 2 7];
fcn_geometry_checkInputsToFunctions(Twocolumn_of_numbers_test, '2column_of_numbers',3);

% Minimum length is 2 or greater
Twocolumn_of_numbers_test = [4 1; 3 9; 2 7];
fcn_geometry_checkInputsToFunctions(Twocolumn_of_numbers_test, '2column_of_numbers',[2 3]);

% Maximum length is 5 or less
Twocolumn_of_numbers_test = [4 1; 3 9; 2 7];
fcn_geometry_checkInputsToFunctions(Twocolumn_of_numbers_test, '2column_of_numbers',[5 4]);

% Maximum length is 3
Twocolumn_of_numbers_test = [4 1; 3 9; 2 7];
fcn_geometry_checkInputsToFunctions(Twocolumn_of_numbers_test, '2column_of_numbers',[3 3]);

assert(true); % pass the test defined by this section if no errors were thrown

%% 2or3column_of_numbers
% 
%   ___           ____            _                               __                       _                   
%  |__ \         |___ \          | |                             / _|                     | |                  
%     ) |___  _ __ __) | ___ ___ | |_   _ _ __ ___  _ __    ___ | |_ _ __  _   _ _ __ ___ | |__   ___ _ __ ___ 
%    / // _ \| '__|__ < / __/ _ \| | | | | '_ ` _ \| '_ \  / _ \|  _| '_ \| | | | '_ ` _ \| '_ \ / _ \ '__/ __|
%   / /| (_) | |  ___) | (_| (_) | | |_| | | | | | | | | || (_) | | | | | | |_| | | | | | | |_) |  __/ |  \__ \
%  |____\___/|_| |____/ \___\___/|_|\__,_|_| |_| |_|_| |_| \___/|_| |_| |_|\__,_|_| |_| |_|_.__/ \___|_|  |___/
%                                                      ______   ______                                         
%                                                     |______| |______|                                        
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Inputs


%% Test the column_of_numbers type (success)
% Test 1 by 2
TwoOrThreeColumn_of_numbers_test = [4 2];
fcn_geometry_checkInputsToFunctions(TwoOrThreeColumn_of_numbers_test, '2or3column_of_numbers');

% Test 1 by 3
TwoOrThreeColumn_of_numbers_test = [4 2 1];
fcn_geometry_checkInputsToFunctions(TwoOrThreeColumn_of_numbers_test, '2or3column_of_numbers');

% Test multiple points - 2 columns
TwoOrThreeColumn_of_numbers_test = [4 1; 3 0; 2 5];
fcn_geometry_checkInputsToFunctions(TwoOrThreeColumn_of_numbers_test, '2or3column_of_numbers');

% Test multiple points - 3 columns
TwoOrThreeColumn_of_numbers_test = [4 1 3; 3 0 5; 2 5 7];
fcn_geometry_checkInputsToFunctions(TwoOrThreeColumn_of_numbers_test, '2or3column_of_numbers');

% Test specified length - 2 columns
TwoOrThreeColumn_of_numbers_test = [4 1; 3 9; 2 7];
fcn_geometry_checkInputsToFunctions(TwoOrThreeColumn_of_numbers_test, '2or3column_of_numbers',3);

% Test specified length - 3 columns
TwoOrThreeColumn_of_numbers_test = [4 1 5; 3 9 5; 2 7 5];
fcn_geometry_checkInputsToFunctions(TwoOrThreeColumn_of_numbers_test, '2or3column_of_numbers',3);

% Minimum length is 2 or greater - 2 columns
TwoOrThreeColumn_of_numbers_test = [4 1; 3 9; 2 7];
fcn_geometry_checkInputsToFunctions(TwoOrThreeColumn_of_numbers_test, '2or3column_of_numbers',[2 3]);

% Minimum length is 2 or greater - 3 columns
TwoOrThreeColumn_of_numbers_test = [4 1 5; 3 9 5; 2 7 5];
fcn_geometry_checkInputsToFunctions(TwoOrThreeColumn_of_numbers_test, '2or3column_of_numbers',[2 3]);

% Maximum length is 5 or less - 2 column
TwoOrThreeColumn_of_numbers_test = [4 1; 3 9; 2 7];
fcn_geometry_checkInputsToFunctions(TwoOrThreeColumn_of_numbers_test, '2or3column_of_numbers',[5 4]);

% Maximum length is 5 or less - 3 column
TwoOrThreeColumn_of_numbers_test = [4 1 4; 3 9 4; 2 7 4];
fcn_geometry_checkInputsToFunctions(TwoOrThreeColumn_of_numbers_test, '2or3column_of_numbers',[5 4]);

% Length MUST be 3 - 2 column
TwoOrThreeColumn_of_numbers_test = [4 1; 3 9; 2 7];
fcn_geometry_checkInputsToFunctions(TwoOrThreeColumn_of_numbers_test, '2or3column_of_numbers',[3 3]);

% Length MUST be 3 - 3 column
TwoOrThreeColumn_of_numbers_test = [4 1 3; 3 9 3; 2 7 3];
fcn_geometry_checkInputsToFunctions(TwoOrThreeColumn_of_numbers_test, '2or3column_of_numbers',[3 3]);

assert(true); % pass the test defined by this section if no errors were thrown

%% Fail conditions
if 1==0
    %% column_of_numbers
    %             _                               __                       _
    %            | |                             / _|                     | |
    %    ___ ___ | |_   _ _ __ ___  _ __    ___ | |_ _ __  _   _ _ __ ___ | |__   ___ _ __ ___
    %   / __/ _ \| | | | | '_ ` _ \| '_ \  / _ \|  _| '_ \| | | | '_ ` _ \| '_ \ / _ \ '__/ __|
    %  | (_| (_) | | |_| | | | | | | | | || (_) | | | | | | |_| | | | | | | |_) |  __/ |  \__ \
    %   \___\___/|_|\__,_|_| |_| |_|_| |_| \___/|_| |_| |_|\__,_|_| |_| |_|_.__/ \___|_|  |___/
    %                                  ______   ______
    %                                 |______| |______|
    %
    %
    
    
    %% Test the column_of_numbers type (FAILURE because 1 x 2)
    column_of_numbers_test = [4 1];
    fcn_geometry_checkInputsToFunctions(column_of_numbers_test, 'column_of_numbers');
    
    assert(true); % pass the test defined by this section if no errors were thrown

    %% Test the column_of_numbers type (FAILURE because 3 long, not 2)
    column_of_numbers_test = [4; 3; 2];
    fcn_geometry_checkInputsToFunctions(column_of_numbers_test, 'column_of_numbers',2);
    
    assert(true); % pass the test defined by this section if no errors were thrown
   
    %% Test the column_of_numbers type (FAILURE because NaN)
    column_of_numbers_test = [4; nan; 2];
    fcn_geometry_checkInputsToFunctions(column_of_numbers_test, 'column_of_numbers',3);
    
    assert(true); % pass the test defined by this section if no errors were thrown

    %% 2column_of_numbers
    %
    %   ___           _                               __                       _
    %  |__ \         | |                             / _|                     | |
    %     ) |___ ___ | |_   _ _ __ ___  _ __    ___ | |_ _ __  _   _ _ __ ___ | |__   ___ _ __ ___
    %    / // __/ _ \| | | | | '_ ` _ \| '_ \  / _ \|  _| '_ \| | | | '_ ` _ \| '_ \ / _ \ '__/ __|
    %   / /| (_| (_) | | |_| | | | | | | | | || (_) | | | | | | |_| | | | | | | |_) |  __/ |  \__ \
    %  |____\___\___/|_|\__,_|_| |_| |_|_| |_| \___/|_| |_| |_|\__,_|_| |_| |_|_.__/ \___|_|  |___/
    %                                      ______   ______
    %                                     |______| |______|
    %
    % See: http://patorjk.com/software/taag/#p=display&f=Big&t=Inputs
    
    %% Test the column_of_numbers type (FAILURE because 1 x 1)
    Twocolumn_of_numbers_test = [4];
    fcn_geometry_checkInputsToFunctions(Twocolumn_of_numbers_test, '2column_of_numbers');
    
    assert(true); % pass the test defined by this section if no errors were thrown

    %% Test the column_of_numbers type (FAILURE because 3 long, not 2)
    Twocolumn_of_numbers_test = [4 1; 3 1; 2 2];
    fcn_geometry_checkInputsToFunctions(Twocolumn_of_numbers_test, '2column_of_numbers',2);
    
    assert(true); % pass the test defined by this section if no errors were thrown
   
    %% Test the column_of_numbers type (FAILURE because NaN)
    Twocolumn_of_numbers_test = [4 1; nan 1; 2 0];
    fcn_geometry_checkInputsToFunctions(Twocolumn_of_numbers_test, '2column_of_numbers',3);
        
    assert(true); % pass the test defined by this section if no errors were thrown

    %% Minimum length is 4 or greater
    Twocolumn_of_numbers_test = [4 1; 3 9; 2 7];
    fcn_geometry_checkInputsToFunctions(Twocolumn_of_numbers_test, '2column_of_numbers',[4 5]);
    
    assert(true); % pass the test defined by this section if no errors were thrown

    %% Maximum length is 2 or less
    Twocolumn_of_numbers_test = [4 1; 3 9; 2 7];
    fcn_geometry_checkInputsToFunctions(Twocolumn_of_numbers_test, '2column_of_numbers',[2 1]);
    
    assert(true); % pass the test defined by this section if no errors were thrown

    %% Maximum length is 2
    Twocolumn_of_numbers_test = [4 1; 3 9; 2 7];
    fcn_geometry_checkInputsToFunctions(Twocolumn_of_numbers_test, '2column_of_numbers',[2 2]);
    
    assert(true); % pass the test defined by this section if no errors were thrown

    %% 2or3column_of_numbers
    %
    %   ___           ____            _                               __                       _
    %  |__ \         |___ \          | |                             / _|                     | |
    %     ) |___  _ __ __) | ___ ___ | |_   _ _ __ ___  _ __    ___ | |_ _ __  _   _ _ __ ___ | |__   ___ _ __ ___
    %    / // _ \| '__|__ < / __/ _ \| | | | | '_ ` _ \| '_ \  / _ \|  _| '_ \| | | | '_ ` _ \| '_ \ / _ \ '__/ __|
    %   / /| (_) | |  ___) | (_| (_) | | |_| | | | | | | | | || (_) | | | | | | |_| | | | | | | |_) |  __/ |  \__ \
    %  |____\___/|_| |____/ \___\___/|_|\__,_|_| |_| |_|_| |_| \___/|_| |_| |_|\__,_|_| |_| |_|_.__/ \___|_|  |___/
    %                                                      ______   ______
    %                                                     |______| |______|
    % See: http://patorjk.com/software/taag/#p=display&f=Big&t=Inputs
    
    

    %% Test the column_of_numbers type (FAILURE because 1 x 1)
    TwoOrThreeColumn_of_numbers_test = [4];
    fcn_geometry_checkInputsToFunctions(TwoOrThreeColumn_of_numbers_test, '2or3column_of_numbers');

    assert(true); % pass the test defined by this section if no errors were thrown

    %% Test the column_of_numbers type (FAILURE because 1 x 4)
    TwoOrThreeColumn_of_numbers_test = [4 1 1 1];
    fcn_geometry_checkInputsToFunctions(TwoOrThreeColumn_of_numbers_test, '2or3column_of_numbers');

    assert(true); % pass the test defined by this section if no errors were thrown

    %% Test the column_of_numbers type (FAILURE because 3 long, not 2)
    TwoOrThreeColumn_of_numbers_test = [4 1; 3 1; 2 2];
    fcn_geometry_checkInputsToFunctions(TwoOrThreeColumn_of_numbers_test, '2or3column_of_numbers',2);
    
    assert(true); % pass the test defined by this section if no errors were thrown
   
    %% Test the column_of_numbers type (FAILURE because NaN)
    TwoOrThreeColumn_of_numbers_test = [4 1; nan 1; 2 0];
    fcn_geometry_checkInputsToFunctions(TwoOrThreeColumn_of_numbers_test, '2or3column_of_numbers',3);
    
    assert(true); % pass the test defined by this section if no errors were thrown

    %% Minimum length is 4 or greater
    TwoOrThreeColumn_of_numbers_test = [4 1; 3 9; 2 7];
    fcn_geometry_checkInputsToFunctions(TwoOrThreeColumn_of_numbers_test, '2or3column_of_numbers',[4 5]);
    
    assert(true); % pass the test defined by this section if no errors were thrown

    %% Maximum length is 2 or less
    TwoOrThreeColumn_of_numbers_test = [4 1; 3 9; 2 7];
    fcn_geometry_checkInputsToFunctions(TwoOrThreeColumn_of_numbers_test, '2or3column_of_numbers',[2 1]);
    
    assert(true); % pass the test defined by this section if no errors were thrown

    %% Maximum length is 2
    TwoOrThreeColumn_of_numbers_test = [4 1; 3 9; 2 7];
    fcn_geometry_checkInputsToFunctions(TwoOrThreeColumn_of_numbers_test, '2or3column_of_numbers',[2 2]);
    
    assert(true); % pass the test defined by this section if no errors were thrown

end
