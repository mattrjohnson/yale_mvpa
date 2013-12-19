function [subs results] = mrj_eeg_classify_doclassification( subs, feature_selection, classifier )

% very simple delegator to either classification over subjects, or k-fold cross-validation, with or without classifying by subject
%  - growing more complex now -- will branch again over whether we have separate training and test sets

if ischar(classifier.kfold) && strcmp(classifier.kfold,'runs') %special value for fMRI Img_Perc
    [subs results] = mrj_eeg_classify_doclassification_kfold_byrun( subs, feature_selection, classifier );
elseif isfield( subs, 'trainset' )
    % for the moment, assume we will only want to be doing training/test within subs instead of between
    %  - also won't bother implementing by_timepoint yet, though we may want to at some point
    if classifier.by_timepoint
        error( 'Not yet implemented!' );
    else
        [subs results] = mrj_eeg_classify_doclassification_traintest( subs, feature_selection, classifier );
    end
elseif classifier.by_timepoint
    if classifier.kfold == 0 % classify over subjects rather than doing k-fold cross-validation within subs
        [subs results] = mrj_eeg_classify_doclassification_over_subjects_bytimepoint( subs, feature_selection, classifier );
    else
        [subs results] = mrj_eeg_classify_doclassification_kfold_bytimepoint( subs, feature_selection, classifier );
    end
else
    if classifier.kfold == 0 % classify over subjects rather than doing k-fold cross-validation within subs
        [subs results] = mrj_eeg_classify_doclassification_over_subjects( subs, feature_selection, classifier );
    else
        [subs results] = mrj_eeg_classify_doclassification_kfold( subs, feature_selection, classifier );
    end
end