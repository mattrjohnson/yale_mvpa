function classifier = yale_mvpa_config_classifier
%------------------------------------------------

classifier.separate_traintest =         0;                  % Choices: 0 or 1
                                                            % Specify 1 if you have separate datasets for training and testing, 0 otherwise. If
                                                            %  this is set to 1, then classifier.kfold (below) should be set to an empty matrix
                                                            %  as the two options are incompatible. If you are using this option, you also need
                                                            %  to make sure that your actual data are split into a training set and a test set.
                                                            %  The built-in EEG pre-processing scripts (when they are finished) will do this if
                                                            %  the correct options are specified, but otherwise you're on your own.

classifier.kfold =                      'runs';             % Choices: Empty, 'runs', 'subs', or a positive integer > 1
                                                            % Should only be empty if specifying separate training and test datasets using the
                                                            %  option above (classifier.separate_traintest). If set to 'runs', then will use a
                                                            %  leave-one-run-out cross-validation approach (mostly for fMRI data, but when the
                                                            %  new EEG pre-processing functions are done, those will include run information
                                                            %  as well. If set to 'subs', will use a leave-one-subject-out cross-validation
                                                            %  approach. If set to an integer (e.g. 5), then will do that fold cross-validation
                                                            %  (e.g. 5-fold). In the integer case, the data will be divvied up randomly, so you
                                                            %  will likely want to run multiple iterations (see classifier.nits below).
                                                            
classifier.nits =                       '';                 % Choices: 0/empty or a positive integer
                                                            % This option only applies if you enter an integer for classifier.kfold above. It
                                                            %  specifies the number of iterations for the classifier to run; the results will be
                                                            %  averaged over all iterations in the end. If not using cross-validation, or if
                                                            %  cross-validating over runs or subjects, this option should be ignored -- but it
                                                            %  would be good practice to set it to 0 or an empty matrix in those cases, just to 
                                                            %  make it more obvious that 'nits' is not being used.
                                                            
                                                            % (or also if separate_traintest and different #s of trials per condition, elaborate
                                                            %  later when that situation is supported)
                                                            
                                                            % EXCEPTION to the general statement above: If you are randomly shuffling your data
                                                            %  to test the classifier (using the classifier.shuffle_data_randomly option... see
                                                            %  directly below), then you can (and probably should) set this value to an integer
                                                            %  larger than 0, so that you can try multiple randomly generated datasets.
                                                            
                                                            % EDIT: Actually, it is not currently true that nits will be ignored, so need to 
                                                            %  rephrase that above. Although we have just made the fix for 0/empty so that is
                                                            %  OK to specify
                                                            
classifier.shuffle_data_randomly =      0;                  % Choices: 0, 1, or 2
                                                            % Should be left at 0 most of the time, which indicates that data will not be 
                                                            %  shuffled; normally this is what you want to do if you want to get actual above-
                                                            %  chance performance. Other values are for testing to make sure the classification
                                                            %  algorithms are not biased; they should yield chance performance if everything is
                                                            %  working correctly. If this option is set to 1, condition labels will be shuffled 
                                                            %  randomly before classification. If it is set to 2, data will be pre-analyzed as 
                                                            %  usual, but when the time to classify comes, the algorithm will just stick in 
                                                            %  completely randomly generated data of the same size.

classifier.cheat =                      0;                  % Choices: 0 or 1
                                                            % Should be left at 0 most of the time, because this option will replace your real
                                                            %  data with fake data that should yield high-accuracy classification, to test your
                                                            %  algorithms. Be careful with this one and if you change it to 1, remember to 
                                                            %  change it back before performing real analyses! (Note: mutually incompatible with
                                                            %  classifier.shuffle_data_randomly above, because this option is basically the
                                                            %  opposite of that.

classifier.cheat_noise_factor =         1;                  % If classifier.cheat is set to 1, this is a multiplier indicating how much noise to
                                                            %  add to the fake data. 0 would indicate no noise, i.e. the data should be perfectly
                                                            %  classifiable with no chance of error (although perhaps some classification 
                                                            %  algorithms might not like overly perfect input data?). 1 would indicate a maximum
                                                            %  noise magnitude equal to the difference between conditions (for a given feature/
                                                            %  trial). Useful values may depend on how much data you're putting in and how noise-
                                                            %  tolerant your classification algorithm is (and this option might be useful for
                                                            %  testing that), but think small, positive numbers, i.e., 1-10ish?

classifier.mean_treatment =             3;                  % Choices: 0, 1, 2, or 3
                                                            % Specifies what to do about the mean signal across features for each trial during
                                                            %  classification. If set to 0, there will be no alteration from the original data.
                                                            %  A setting of 1 specifies that the mean across features will be subtracted out
                                                            %  (i.e., each trial's mean will become 0). A setting of 2 specifies that ONLY the
                                                            %  mean will be kept -- thus each trial will only have a single feature, which is
                                                            %  the mean of all the original features. A setting of 3 specifies that each trial
                                                            %  will be z-scored across features before classification -- meaning that not only
                                                            %  the mean of all features (on a given trial) will be 0, but the standard deviation
                                                            %  will also be 1.


% which classifier to use; classifier.args will depend on what your train and test functions are
    % Matlab SVM
classifier.trainfunc =                  @yale_mvpa_train_matlabsvm;
classifier.testfunc =                   @yale_mvpa_test_matlabsvm;
classifier.args =                       [];
    % Logistic regression
% classifier.trainfunc =                  @train_logreg;
% classifier.testfunc =                   @test_logreg;
% classifier.args.penalty =               100;
    % SMLR
% classifier.trainfunc =                  @train_smlr;
% classifier.testfunc =                   @test_smlr;
% classifier.args.verbose =               0;
    % Ridge regression
% classifier.trainfunc =                  @train_ridge;
% classifier.testfunc =                   @test_ridge;
% classifier.args.penalty =               100;





% add in an option to override basic classification functions (e.g. if you want to write your own thing, like a pattern similarity analyzer)
% note about timepoint flag in general settings once that's added back in