function [subs results] = mrj_eeg_classify_doclassification_kfold_bytimepoint( subs, feature_selection, classifier )

n_subs =                                length(subs);
results =                               cell( n_subs, 1 );
if classifier.do_parallel
    parfor i=1:n_subs
        results{i} =                        mrj_eeg_classify_doclassification_kfold_bytimepoint_onesub( subs(i), feature_selection, classifier, i );
    end
else
    for i=1:n_subs
        results{i} =                        mrj_eeg_classify_doclassification_kfold_bytimepoint_onesub( subs(i), feature_selection, classifier, i );
    end
end

%--------------------------------------------------------------------------------------------------------------------------
%--------------------------------------------------------------------------------------------------------------------------
function results = mrj_eeg_classify_doclassification_kfold_bytimepoint_onesub( sub, feature_selection, classifier, subnum )

% sub{i}.classifer_data should still be a cell, but each element of the cell should also be a cell of timepoints
n_timepoints =                          length( sub.classifier_data{1} );

% initialize output struct
% n.b.: now results is going to be a struct arry with n_timepoints elements, rather than a single-element struct
results(n_timepoints).n_features =      zeros( classifier.nits, 1 );
results(n_timepoints).acts =            cell( classifier.nits, 1 );
results(n_timepoints).testtargs =       cell( classifier.nits, 1 );

% figure out some preliminary classifier stuff that is common across timepoints
n_trials_min =                          min(sub.n_trials_per_cond);
n_testtrials_percond =                  floor(n_trials_min/classifier.kfold);
n_traintrials_percond =                 n_testtrials_percond*(classifier.kfold-1);
% traintargs and testtargs are condition labels for train and test, respectively
traintargs =                            kron( eye(sub.n_conds), ones(1,n_traintrials_percond) );
testtargs =                             kron( eye(sub.n_conds), ones(1,n_testtrials_percond) );
sub_tmp.n_conds =                       sub.n_conds;

% loop over timepoints, then iterations; different timepoints may have different numbers of features, etc.
for i=1:n_timepoints % loop over timepoints
    % display update
    disp(['Subject ' int2str(subnum) ', timepoint ' int2str(i) ' of ' int2str(n_timepoints)]);
    
    for j=1:sub_tmp.n_conds % loop over conditions, set up classifier info for this timepoint
        sub_tmp.classifier_data{j} =    sub.classifier_data{j}{i};
    end
    n_features =                        size( sub_tmp.classifier_data{1}, 1 );
    
    % do actual classification
    for j=1:classifier.nits
        % shuffle condition labels?
        if classifier.shuffle_data_randomly==1 % shuffles condition labels randomly
            traintargs =                traintargs(:,randperm(size(traintargs,2)));
            testtargs =                 testtargs(:,randperm(size(testtargs,2)));
        end
        
        % divide classifier data into training and test sets
        %   n.b.: if classifier.shuffle_data_randomly is set to 2, that happens in this function too
        [trainpats testpats] = mrj_eeg_classify_doclassification_kfold_divide_traintest( sub_tmp, classifier, n_features, n_traintrials_percond, n_testtrials_percond, traintargs, testtargs );
    
        % do feature selection, if specified
        if feature_selection.use
            feature_inds =              feval( feature_selection.function, trainpats, traintargs, feature_selection.args );
            trainpats =                 trainpats( feature_inds, : );
            testpats =                  testpats( feature_inds,: );
            disp([' - feature selection: ' int2str(size(trainpats,1)) ' features']);
        end
        results(i).n_features(j) =      size(trainpats,1);
        
        % do actual classification here
        s =                             feval( classifier.trainfunc, trainpats, traintargs, classifier.args );
        [acts s] =                      feval( classifier.testfunc, testpats, testtargs, s );
        
        % any special behavior for particular classification functions can go here
        if isequal(classifier.testfunc,@test_matlabsvm_mrj)
            results(i).svm_orig_accs{j} =   [ s(:).accs ];
        end
        
        results(i).acts{j} =            acts;
        results(i).testtargs{j} =       testtargs; % testtargs will generally be the same for all iterations, but not necessarily
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% KEY SUB-FUNCTION(S)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
function [trainpats testpats] = mrj_eeg_classify_doclassification_kfold_divide_traintest( sub, classifier, n_features, n_traintrials_percond, n_testtrials_percond, traintargs, testtargs )
% note: if I did this right, this function should need no modifications from regular version to by_timepoint version, other than this comment

trainpats =                         zeros( n_features, n_traintrials_percond * sub.n_conds );
testpats =                          zeros( n_features, n_testtrials_percond * sub.n_conds );

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
    end
end
