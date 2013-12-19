function [subs results] = mrj_eeg_classify_doclassification_kfold( subs, feature_selection, classifier )

n_subs =                                length(subs);
results =                               cell( n_subs, 1 );
if classifier.do_parallel
    parfor i=1:n_subs
        results{i} =                    mrj_eeg_classify_doclassification_kfold_onesub( subs(i), feature_selection, classifier, i );
    end
else
    for i=1:n_subs
        results{i} =                    mrj_eeg_classify_doclassification_kfold_onesub( subs(i), feature_selection, classifier, i );
    end
end

%--------------------------------------------------------------------------------------------------------------
%--------------------------------------------------------------------------------------------------------------
function results = mrj_eeg_classify_doclassification_kfold_onesub( sub, feature_selection, classifier, subnum )

%added for Img_Perc fMRI
sub.n_trials_per_cond = size(sub.classifier_data{1},2);

% initialize output struct
results.n_features =                    zeros( classifier.nits, 1 );
results.acts =                          cell( classifier.nits, 1 );
results.testtargs =                     cell( classifier.nits, 1 );
results.traininds =                     cell( classifier.nits, 1 );
results.testinds =                      cell( classifier.nits, 1 );
results.train_eegfile_trial_inds =      cell( classifier.nits, 1 );
results.test_eegfile_trial_inds =       cell( classifier.nits, 1 );

% figure out some preliminary classifier stuff
n_features =                            size( sub.classifier_data{1}, 1 );
n_trials_min =                          min(sub.n_trials_per_cond);
n_testtrials_percond =                  floor(n_trials_min/classifier.kfold);
n_traintrials_percond =                 n_testtrials_percond*(classifier.kfold-1);
% traintargs and testtargs are condition labels for train and test, respectively
traintargs =                            kron( eye(sub.n_conds), ones(1,n_traintrials_percond) );
testtargs =                             kron( eye(sub.n_conds), ones(1,n_testtrials_percond) );

% do actual classification
for i=1:classifier.nits
    % display update
    disp(['Subject ' int2str(subnum) ', iteration ' int2str(i) ' of ' int2str(classifier.nits)]);
    
    % shuffle condition labels?
    if classifier.shuffle_data_randomly==1 % shuffles condition labels randomly
        traintargs =                    traintargs(:,randperm(size(traintargs,2)));
        testtargs =                     testtargs(:,randperm(size(testtargs,2)));
    end
    
    % divide classifier data into training and test sets
    %   n.b.: if classifier.shuffle_data_randomly is set to 2, that happens in this function too
    [trainpats testpats traininds testinds] = mrj_eeg_classify_doclassification_kfold_divide_traintest( sub, classifier, n_features, n_traintrials_percond, n_testtrials_percond, traintargs, testtargs );

    % do feature selection, if specified
    if feature_selection.use
        feature_inds =                  feval( feature_selection.function, trainpats, traintargs, feature_selection.args );
        trainpats =                     trainpats( feature_inds, : );
        testpats =                      testpats( feature_inds,: );
        disp([' - feature selection: ' int2str(size(trainpats,1)) ' features']);
    end
    results.n_features(i) =             size(trainpats,1);
    
    % do actual classification here
    s =                                 feval( classifier.trainfunc, trainpats, traintargs, classifier.args );
    [acts s] =                          feval( classifier.testfunc, testpats, testtargs, s );
    
    % any special behavior for particular classification functions can go here
    if isequal(classifier.testfunc,@test_matlabsvm_mrj)
        results.svm_orig_accs{i} =      [ s(:).accs ];
    end
    
    results.acts{i} =                   acts;
    results.testtargs{i} =              testtargs; % testtargs will generally be the same for all iterations, but not necessarily
    results.traininds{i} =              traininds;
    results.testinds{i} =               testinds;
    
    %commented out for Img_Perc fMRI
%     train_eegfile_trial_inds_tmp =      zeros( size(traininds,2), 2 );
%     test_eegfile_trial_inds_tmp =       zeros( size(testinds,2),  2 );
%     for j=1:size( traininds, 1 ) %loop thru conditions
%         traininds_tmp =                 traininds(j,:)>0;
%         testinds_tmp =                  testinds(j,:)>0;
%         
%         train_eegfile_trial_inds_tmp( traininds_tmp, 1 ) = ...
%                                         sub.eegfile_inds{j}( traininds(j,traininds_tmp) );
%         train_eegfile_trial_inds_tmp( traininds_tmp, 2 ) = ...
%                                         sub.trial_inds{j}( traininds(j,traininds_tmp) );
%         test_eegfile_trial_inds_tmp(  testinds_tmp, 1 )  = ...
%                                         sub.eegfile_inds{j}( testinds(j,testinds_tmp)   );
%         test_eegfile_trial_inds_tmp(  testinds_tmp, 2 )  = ...
%                                         sub.trial_inds{j}( testinds(j,testinds_tmp)   );
%     end
%     results.train_eegfile_trial_inds{i} = ...
%                                         train_eegfile_trial_inds_tmp;
%     results.test_eegfile_trial_inds{i} = ...
%                                         test_eegfile_trial_inds_tmp;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% KEY SUB-FUNCTION(S)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
function [trainpats testpats all_traininds all_testinds] = mrj_eeg_classify_doclassification_kfold_divide_traintest( sub, classifier, n_features, n_traintrials_percond, n_testtrials_percond, traintargs, testtargs )

trainpats =                         zeros( n_features, n_traintrials_percond * sub.n_conds );
testpats =                          zeros( n_features, n_testtrials_percond * sub.n_conds );
all_traininds =                     zeros( size(traintargs) );
all_testinds =                      zeros( size(testtargs) );

if classifier.shuffle_data_randomly == 2 %just stick in random data the same size as the real data
    trainpats =                     rand(size(trainpats));
    testpats =                      rand(size(testpats));
else
    % loop through conditions; divide data into training and test sets
    for i=1:sub.n_conds
        inds =                          randperm(size(sub.classifier_data{i},2));
        traininds =                     inds(1:n_traintrials_percond);
        testinds =                      inds((n_traintrials_percond+1):(n_traintrials_percond+n_testtrials_percond));
        trainpats(:,logical(traintargs(i,:))) = ...
                                        sub.classifier_data{i}(:,traininds); %indexed nfeatures x ntrials
        testpats(:,logical(testtargs(i,:))) = ...
                                        sub.classifier_data{i}(:,testinds); %indexed nfeatures x ntrials
        all_traininds(i, logical(traintargs(i,:))) = ...
                                        traininds;
        all_testinds(i, logical(testtargs(i,:))) = ...
                                        testinds;
    end
end
