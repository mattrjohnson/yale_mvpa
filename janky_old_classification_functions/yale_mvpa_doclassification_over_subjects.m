function [subs results] = mrj_eeg_classify_doclassification_over_subjects( subs, feature_selection, classifier )

n_subs =                                length(subs);
results =                               cell( n_subs, 1 );
if classifier.do_parallel
    parfor i=1:n_subs
        results{i} =                        mrj_eeg_classify_doclassification_over_subjects_onesub( subs, feature_selection, classifier, i );
    end
else
    for i=1:n_subs
        results{i} =                        mrj_eeg_classify_doclassification_over_subjects_onesub( subs, feature_selection, classifier, i );
    end
end

%-------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------
function results = mrj_eeg_classify_doclassification_over_subjects_onesub( subs, feature_selection, classifier, test_sub )

%initialize output struct
results.n_features =                    zeros( classifier.nits, 1 );
results.acts =                          cell( classifier.nits, 1 );
results.testtargs =                     cell( classifier.nits, 1 );

% figure out some preliminary classifier stuff
n_subs =                                length(subs);
n_trials_persub_percond =               [subs.n_trials_per_cond]'; % indexed by subject, then condition
training_subs =                         setdiff( 1:n_subs, test_sub );
n_traintrials_min =                     min(min( n_trials_persub_percond( training_subs, : ) ));
n_testtrials_min =                      min( n_trials_persub_percond( test_sub, : ) );
n_features =                            size( subs(1).classifier_data{1}, 1 ); % assume all subs have same # of features?
traintargs =                            kron( eye(subs(1).n_conds), ones(1,n_traintrials_min * (n_subs-1)) );
testtargs =                             kron( eye(subs(1).n_conds), ones(1,n_testtrials_min) );

% do actual classification
for i=1:classifier.nits
    % display update
    disp(['Subject ' int2str(test_sub) ', iteration ' int2str(i) ' of ' int2str(classifier.nits)]);
    
    % shuffle condition labels?
    if classifier.shuffle_data_randomly==1 % shuffles condition labels randomly
        traintargs =                    traintargs(:,randperm(size(traintargs,2)));
        testtargs =                     testtargs(:,randperm(size(testtargs,2)));
    end

    % divide classifier data into training and test sets
    %   n.b.: if classifier.shuffle_data_randomly is set to 2, that happens in this function too
    [trainpats testpats] = mrj__doclassification_over_subjects_divide_traintest( subs, classifier, n_features, n_traintrials_min, n_testtrials_min, testtargs, test_sub );
    
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
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% KEY SUB-FUNCTION(S)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
function [trainpats testpats] = mrj__doclassification_over_subjects_divide_traintest( subs, classifier, n_features, n_traintrials_percond, n_testtrials_percond, testtargs, test_sub )
% shortenened name due to length constraints

n_subs =                            length( subs );
trainpats =                         zeros( n_features, n_traintrials_percond * subs(1).n_conds * (n_subs - 1) );
testpats =                          zeros( n_features, n_testtrials_percond * subs(1).n_conds );
training_subs =                     setdiff( 1:n_subs, test_sub );

if classifier.shuffle_data_randomly == 2 %just stick in random data the same size as the real data
    trainpats =                     rand(size(trainpats));
    testpats =                      rand(size(testpats));
else
    % loop through conditions; divide data into training and test sets
    for i=1:subs(1).n_conds
        % first do test data for the test subjet
        inds =                          randperm(size(subs(test_sub).classifier_data{i},2));
        inds =                          inds( 1:n_testtrials_percond );
        testpats(:,logical(testtargs(i,:))) = ...
                                        subs(test_sub).classifier_data{i}(:,inds);
        
        % now loop through training subjects
        for j=1:length(training_subs)
            this_sub =                  training_subs(j);
            inds =                      randperm(size(subs(this_sub).classifier_data{i},2));
            inds =                      inds( 1:n_traintrials_percond );
            start_ind =                 (i-1)*(n_subs-1)*(n_traintrials_percond) + (j-1)*n_traintrials_percond + 1;
            end_ind =                   start_ind + n_traintrials_percond - 1;
            trainpats(:,start_ind:end_ind) = ...
                                        subs(this_sub).classifier_data{i}(:,inds);
        end
    end
end
