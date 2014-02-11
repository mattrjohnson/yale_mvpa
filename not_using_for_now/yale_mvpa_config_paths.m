function paths = yale_mvpa_config_paths
%--------------------------------------


%Took the general path function out because that was the only thing in here getting used. It's now in _config_general.
% For EEG-processing-specific path stuff, it makes more sense to put that in the EEG-processing plugin settings
% So when we re-implement EEG processing, we can stick these back in.




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
