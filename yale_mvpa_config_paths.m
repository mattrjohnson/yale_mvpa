function paths = yale_mvpa_config_paths
%--------------------------------------


% paths.path_function =                   'mrj_eeg_paths_memlabs';
% paths.path_function = 'addpath /Users/matt/Desktop/docs/data/Img_Perc/from_refnegprime_tmp_scripts';
paths.path_function =                   '';                 % String describing a function to be eval'd in order to add any appropriate paths
                                                            %  (or do any other initialization/setup). Can be anything you'd enter on the Matlab
                                                            %  prompt. The toolbox will take care of its own paths but you may want this for other
                                                            %  home-grown functions and so forth. Leave as an empty string if you don't need it.



% all this stuff is theoretically specific to EEG data, let's see how much we can get away with removing for now

% paths.experiment_directory =            '';
% paths.p_arf =                           '';
% paths.p_mon =                           '';
% 
% paths.sub_list =                        {  };
% % paths.eeg_subdir =                      '';
% % paths.p_sgc =                           '';
% 
%     % training and testing on different data
% paths.train.eeg_subdir =                '';
% paths.train.p_sgc =                     { };
% paths.test.eeg_subdir =                 '';
% paths.test.p_sgc =                      { };
