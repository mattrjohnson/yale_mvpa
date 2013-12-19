function [subs results] = mrj_eeg_classify_doclassification_kfold_byrun( subs, feature_selection, classifier )
% special function for fMRI Img_Perc

n_subs =                                length(subs);
results =                               cell( n_subs, 1 );
if classifier.do_parallel
    parfor i=1:n_subs
        results{i} =                    mrj_eeg_classify_doclassification_kfold_byrun_onesub( subs(i), feature_selection, classifier, i );
    end
else
    for i=1:n_subs
        results{i} =                    mrj_eeg_classify_doclassification_kfold_byrun_onesub( subs(i), feature_selection, classifier, i );
    end
end

%--------------------------------------------------------------------------------------------------------------
%--------------------------------------------------------------------------------------------------------------
function results = mrj_eeg_classify_doclassification_kfold_byrun_onesub( sub, feature_selection, classifier, subnum )

%normally, sub.classifier_data is a cell vector (by condition), each containing an nfeatures x ntrials matrix
% for Img_Perc, sub.classifier_data will be a cell matrix, n_conds x n_runs

if ~classifier.shuffle_data_randomly
    classifier.nits = 1;
end
n_runs = size(sub.classifier_data,2);
n_trials_percond_perrun = size(sub.classifier_data{1},2);
n_conds = sub.n_conds;

% initialize output struct
results.n_features =                    zeros( classifier.nits, 1 );
results.acts =                          cell( classifier.nits, 1 );
results.testtargs =                     cell( classifier.nits, 1 );
% results.traininds =                     cell( classifier.nits, 1 ); %won't actually get used
% results.testinds =                      cell( classifier.nits, 1 ); %won't actually get used
% results.train_eegfile_trial_inds =      cell( classifier.nits, 1 ); %won't actually get used
% results.test_eegfile_trial_inds =       cell( classifier.nits, 1 ); %won't actually get used

% figure out some preliminary classifier stuff
n_features =                            size( sub.classifier_data{1}, 1 );
% n_trials_min =                          min(sub.n_trials_per_cond);
n_testtrials_percond =                  n_trials_percond_perrun; %changed for Img_Perc
% n_traintrials_percond =                 n_trials_percond_perrun * (n_runs - 1); %changed for Img_Perc %don't really need anymore, actually
% traintargs and testtargs are condition labels for train and test, respectively
%traintargs =                            kron( eye(sub.n_conds), ones(1,n_traintrials_percond) ); %changing for Img_Perc
testtargs =                             kron( eye(sub.n_conds), ones(1,n_testtrials_percond) );
for i=1:(n_runs-1)
    traintargs_3d(:,:,i) =              kron( eye(sub.n_conds), ones(1,n_trials_percond_perrun)); %#ok<AGROW>
                                        %traintargs_3d is n_conds x n_trials_perrun x n_runs
end
traintargs =                            reshape( traintargs_3d, [size(traintargs_3d,1) size(traintargs_3d,2)*size(traintargs_3d,3)] );

%traintargs will be arranged a bit differently for Img_Perc; will contain sections by run instead of having each condition contiguous
% e.g., instead of 
%[1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
% 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
% 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0;
% 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1]

% it'll be 

%[1 1 0 0 0 0 0 0 1 1 0 0 0 0 0 0 1 1 0 0 0 0 0 0 1 1 0 0 0 0 0 0;
% 0 0 1 1 0 0 0 0 0 0 1 1 0 0 0 0 0 0 1 1 0 0 0 0 0 0 1 1 0 0 0 0;
% 0 0 0 0 1 1 0 0 0 0 0 0 1 1 0 0 0 0 0 0 1 1 0 0 0 0 0 0 1 1 0 0;
% 0 0 0 0 0 0 1 1 0 0 0 0 0 0 1 1 0 0 0 0 0 0 1 1 0 0 0 0 0 0 1 1]

%added for Img_Perc
classifier_data_reshaped =              nan(n_conds, n_runs, n_features, n_trials_percond_perrun);
for i=1:size(sub.classifier_data,1) %loop thru conds
    for j=1:size(sub.classifier_data,2) %loop thru runs
        classifier_data_reshaped(i,j,:,:) =                 sub.classifier_data{i,j};
    end
end
classifier_data_reshaped =              permute(classifier_data_reshaped,[2 3 4 1]); %now runs x features x trial_within_cond x cond
classifier_data_reshaped =              reshape(classifier_data_reshaped,[n_runs n_features n_conds*n_trials_percond_perrun]); %now runs x features x trial_within_run

% do actual classification
for i=1:classifier.nits
    % display update
    disp(['Subject ' int2str(subnum) ', iteration ' int2str(i) ' of ' int2str(classifier.nits)]);
    
    % shuffle condition labels? %this whole chunk changed for Img_Perc
    if classifier.shuffle_data_randomly==1 % shuffles condition labels randomly
        for j = 1:n_runs
            traintargs_3d(:,:,j) =      traintargs_3d(:,randperm(size(traintargs_3d,2)),j); %#ok<AGROW>
        end
        traintargs =                    reshape( traintargs_3d, [size(traintargs_3d,1) size(traintargs_3d,2)*size(traintargs_3d,3)] );
        testtargs =                     testtargs(:,randperm(size(testtargs,2)));
    end
    if classifier.shuffle_data_randomly==2
        error('Not yet implemented');
    end
    
    % divide classifier data into training and test sets
    %   n.b.: if classifier.shuffle_data_randomly is set to 2, that happens in this function too
    %[trainpats testpats traininds testinds] = mrj_eeg_classify_doclassification_kfold_divide_traintest( sub, classifier, n_features, n_traintrials_percond, n_testtrials_percond, traintargs, testtargs );
    
    % for Img_Perc, instead of the above, we now need to loop through runs
    acts =                              nan( n_conds, n_conds*n_trials_percond_perrun, n_runs );
    for j=1:n_runs
        train_runs =                    setdiff( 1:n_runs, j );
        trainpats =                     nan( n_features, n_conds*n_trials_percond_perrun, n_runs-1 ); %indexed features x trials_perrun x run
        for k=1:length(train_runs)
            trainpats(:,:,k) =          squeeze( classifier_data_reshaped(train_runs(k),:,:) );
        end
        trainpats =                     reshape( trainpats, [n_features n_conds*n_trials_percond_perrun*(n_runs-1)] );
        testpats =                      squeeze( classifier_data_reshaped( j,:,: ) );
        
        % probably won't need to do feature selection for Img_Perc, so don't worry about it right now
        % do feature selection, if specified
    %     if feature_selection.use
    %         feature_inds =                  feval( feature_selection.function, trainpats, traintargs, feature_selection.args );
    %         trainpats =                     trainpats( feature_inds, : );
    %         testpats =                      testpats( feature_inds,: );
    %         disp([' - feature selection: ' int2str(size(trainpats,1)) ' features']);
    %     end
        results.n_features(i) =             size(trainpats,1); %will need to fix this if we do actually do feature selection
        
        % do actual classification here
        s =                                 feval( classifier.trainfunc, trainpats, traintargs, classifier.args );
        [acts(:,:,j) s] =                   feval( classifier.testfunc, testpats, testtargs, s );
        
        % any special behavior for particular classification functions can go here
        if isequal(classifier.testfunc,@test_matlabsvm_mrj)
            svm_orig_accs(:,:,j) =          [ s(:).accs ]; %#ok<AGROW> %will need to check this to make sure we got ndims / indexing right
        end
    end
    
%         % any special behavior for particular classification functions can go here
%         if isequal(classifier.testfunc,@test_matlabsvm_mrj)
%             results.svm_orig_accs{i} =      [ s(:).accs ];
%         end
    %once we figure out above, will also have to do something here to save results.svm_orig_accs
    
    acts =                              reshape( acts, [n_conds n_conds*n_trials_percond_perrun*n_runs] ); %new line for Img_Perc; now n_conds x n_trials (total, across all runs) again
    
    results.acts{i} =                   acts;
    results.testtargs{i} =              repmat(testtargs,[1 n_runs]); % testtargs will generally be the same for all iterations, but not necessarily
    
    %don't think we need any of the below for Img_Perc
%     results.traininds{i} =              traininds;
%     results.testinds{i} =               testinds;
%     
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


