function results = yale_mvpa_classify( yale_mvpa_config )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONFIGURATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<1 || isempty( yale_mvpa_config )
    yale_mvpa_config.general =           	yale_mvpa_config_general;
    yale_mvpa_config.files =              	yale_mvpa_config_files;
    % yale_mvpa_config.erps =              	mrj_eeg_classify_config_erps;
    % yale_mvpa_config.freqs =               	mrj_eeg_classify_config_freqs;
    % yale_mvpa_config.feature_selection =   	yale_mvpa_config_feature_selection;
    yale_mvpa_config.classifier =        	yale_mvpa_config_classifier;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OTHER SETUP/BOOKKEEPING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% keep track of how long computations take
yale_mvpa_config.start_time.str =     	datestr(now);
yale_mvpa_config.start_time.tictoc =   	tic;

% save original working directory
yale_mvpa_config.owd =                  pwd;

% set paths
yale_mvpa_toolbox_path = fileparts(mfilename('fullpath'));
disp('Adding directories to path (just for the current Matlab session):');
disp(['- ' yale_mvpa_toolbox_path]);
addpath(yale_mvpa_toolbox_path);
disp(['- ' fullfile(yale_mvpa_toolbox_path,'helpers')]);
addpath(fullfile(yale_mvpa_toolbox_path,'helpers'));
if ~isempty( yale_mvpa_config.general.path_function )
    try
        eval( yale_mvpa_config.general.path_function );
    catch %#ok<CTCH>
        error( 'Could not evaluate general.path_function!' );
    end
end

% open parallel computing?
if yale_mvpa_config.general.try_parallel
    try
        disp('Trying to start Matlab pool...');
        matlabpool open;
    catch myexception
        if strcmp(myexception.identifier,'distcomp:interactive:OpenConnection') % pool already open
            disp('Matlab pool appears to be open already, continuing');
        else
            disp('Could not open Matlab pool -- continuing in single-processor mode');
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GET DATA -- VIA EITHER PRE-PROCESSING OR LOADING IN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[yale_mvpa_config, subs] = yale_mvpa_get_data( yale_mvpa_config );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ACTUAL CLASSIFICATION (or other not-technically-classification MVPA)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[yale_mvpa_config, subs, results] = yale_mvpa_doclassification( yale_mvpa_config, subs );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RESULTS (DISPLAY AND/OR SAVE) AND CLEANUP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd(yale_mvpa_config.owd);

%need to modify classification functions to define results functions
% all valid results-processing (and displaying) functions should take 3 parameters and return 1
results = feval(yale_mvpa_config.results.function, yale_mvpa_config, subs, results );

if yale_mvpa_config.files.save_results
    do_save = yale_mvpa_check_save_ok( yale_mvpa_config.files.save_results_fname );
    if do_save
        save( yale_mvpa_config.files.save_results_fname, 'yale_mvpa_config', 'results' );
    end
end

disp('yale_mvpa_classify finished');
toc(yale_mvpa_config.start_time.tictoc);



%later on: maybe add the ability to pass in a directory with configuration files, put the default configuration files in their own
% sub-directory, and only add it to the path if the user does not pass in a directory with their own config files. (although maybe 
% this has been obviated by the change of allowing one to pass in the yale_mvpa_config structure)
