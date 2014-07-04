function files = yale_mvpa_config_files
%--------------------------------------


files.save_preprocessed_data =          0;                  % Choices: 0 or 1
                                                            % If set to 1, will try to save the results of pre-procesing in a file you specify
                                                            %  (via the files.save_preprocessed_data_fname option below). This currently only
                                                            %  applies to EEG data -- fMRI data must be loaded in from file. It also implies
                                                            %  that if this option is set to 1, you must turn on some kind of EEG pre-processing
                                                            %  (ERPs or time-frequency or both). If set to 0, then the option below will be
                                                            %  ignored, but it is probably good practice to set it to an empty string to reduce
                                                            %  confusion.

files.save_preprocessed_data_fname =    '';                 
                                                            % String specifying the filename to write. Can be either an absolute filename, or
                                                            %  relative to the directory that was the Matlab working directory when the user
                                                            %  ran the yale_mvpa_classify script. Must put something in here if the above option
                                                            %  (files.save_preprocessed_data) is set to 1. Leave empty if the option above is
                                                            %  set to 0.
                                                            
files.load_preprocessed_data =          1;                  % Choices: 0 or 1
                                                            % If set to 1, will load in pre-processed data from a file. Note that this will skip
                                                            %  all pre-processing steps for EEG data, even if the settings for those steps would
                                                            %  seem to indicate that they will be run. It will also replace the configuration
                                                            %  options specified for those steps with ones loaded from file, specifying how the
                                                            %  loaded-in data was originally pre-processed. If set to 0, nothing will be loaded,
                                                            %  so you must turn on some kind of EEG pre-processing (ERPs or time-frequency or
                                                            %  both). For fMRI data, currently you MUST specify 1 -- there are currently no 
                                                            %  functions for fMRI pre-processing in this toolbox, so fMRI data must be loaded 
                                                            %  from file.
                                                            % Note that this option is incompatible with the files.save_preprocessed_data 
                                                            %  option... both of them cannot be set to 1 at the same time (although it is OK if
                                                            %  they are both set to 0). If this option is set to 1, then you MUST specify a 
                                                            %  filename in the files.load_preprocessed_data_fname option below. If this option 
                                                            %  is set to 0, then the value in files.load_preprocessed_data_fname will be
                                                            %  ignored, but it is probably good practice to set it to an empty string to reduce
                                                            %  confusion.

files.load_preprocessed_data_fname =    '/Users/matt/Desktop/docs/data/Img_Perc/from_lab_pc/converted_classifier_data/pvsi/preprocessed_data_ImgPercAud_ffa_80voxeach_pvsi.mat';
                                                            % String specifying the filename from which data will be loaded. Can be either an 
                                                            %  absolute filename, or relative to the directory that was the Matlab working 
                                                            %  directory when the user ran the yale_mvpa_classify script. Must put something in 
                                                            %  here if the above option (files.load_preprocessed_data) is set to 1. Leave empty 
                                                            %  if the option above is set to 0.

files.load_data_via_script =            0;                  % Choices: 0 or 1
                                                            % If set to 1, will run a script you specify. (Update more later. Including update
                                                            %  descriptions for above options that are mutually incompatible with this.) OK to
                                                            %  co-occur with save_results parameter below; e.g. if your data-loading function is
                                                            %  some kind of pre-processing function. (In a later update: Change so that the Yale
                                                            %  EEG preprocessing scripts fold into this, and make either load_data_via_script or
                                                            %  load_preprocessed_data mandatory.)

files.load_data_via_script_mfile =      '';                 
                                                            % String specifying the Matlab M-file to run if files.load_data_via_script is set to
                                                            %  1. Should be the name of an M-file (without the actual .m extension) that is 
                                                            %  located somewhere on the Matlab path (including the yale_mvpa directory and any 
                                                            %  paths added by the path function in yale_mvpa_config_general.m). The M-file 
                                                            %  specified should include a function that returns one output, which should match 
                                                            %  the 'subs' data structure used everywhere else (add more detail later).

files.save_results =                    0;                  % Choices: 0 or 1
                                                            % If set to 1, will save the results of MVPA in a file you specify using the 
                                                            %  files.save_results_fname option below. The exact format will depend on the type
                                                            %  of analysis you run. Results should be displayed to the Matlab command window
                                                            %  whether you save them or not. Pretty much the same deal as the options above...
                                                            %  if set to 1, you need to specify a filename in the option below; if set to 0,
                                                            %  good practice to leave the below option empty; etc., etc.

files.save_results_fname =              '';                 
                                                            % String specifying the filename in which to (optionally) write out results. Pretty
                                                            %  much the same deal as the other filename options above (absolute or relative, 
                                                            %  must have something in here if above option is set to 1, leave empty if above
                                                            %  option is 0, etc.).