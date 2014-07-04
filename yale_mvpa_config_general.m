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

general.path_function =                 '';                 % String describing a function to be eval'd in order to add any appropriate paths
                                                            %  (or do any other initialization/setup). Can be anything you'd enter on the Matlab
                                                            %  prompt. The toolbox will take care of its own paths but you may want this for other
                                                            %  home-grown functions and so forth. Leave as an empty string if you don't need it.

general.results_verbose =               1;                  % Choices: 0 or 1
                                                            % If set to 1 (the default), will print out results to the Matlab command window. If
                                                            %  set to 0, results aren't printed at all, so hopefully you're saving them somehow
                                                            %  (either to a file or catching the return value of yale_mvpa_classify)... or else
                                                            %  there isn't much point to having run the classification, is there?


% Still to come: Option to do classification separately by "channel" (fMRI: == ROI?)
% - this option will be mutually exclusive with by_timepoint -- they can't both be used