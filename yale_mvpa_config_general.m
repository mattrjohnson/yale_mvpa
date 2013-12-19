function general = yale_mvpa_config_general
%------------------------------------------

% A few general settings (General Settings! [salute])

general.modality =                      'fmri';             % Choices: 'eeg' or 'fmri'
                                                            % In most places, settings should be modality-independent, but there will be some
                                                            %  situations where it will be a good idea to check and confirm that what we're 
                                                            %  doing makes sense for the given modality.

general.by_timepoint =                  0;                  % Choices: 0 or 1
                                                            % Single universal flag indicating whether pre-processing and classification will 
                                                            %  proceed over timepoints, rather than omnibus classification. Originally intended
                                                            %  for EEG but can work for fMRI as well.

general.try_parallel =                  0;                  % Choices: 0 or 1
                                                            % If set to 1, will try to do processing in parallel (using Matlab Parallel Computing
                                                            %  Toolbox) wherever possible.

% Still to come: Option to do classification separately by "channel" (fMRI: == ROI?)
% - this option will be mutually exclusive with by_timepoint -- they can't both be used