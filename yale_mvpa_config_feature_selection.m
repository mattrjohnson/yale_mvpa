function feature_selection = mrj_eeg_classify_config_feature_selection
%---------------------------------------------------------------------


feature_selection.use =                 0;
% which feature selection to use; feature_selection.args will depend on what your feature_selection.function is
    % simple ANOVA over training data
feature_selection.function =            @mrj_feature_selection_anova;
feature_selection.args.thresh =         .01;
feature_selection.args.parallel =       0;                  % only set this == 1 if classifier.do_parallel is 0
                                                            % if set, will use a parfor loop in feature selection even if not running main classification in parallel
