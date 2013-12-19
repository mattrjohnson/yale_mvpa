function [subs results] = mrj_eeg_classify_doclassification_traintest( subs, feature_selection, classifier )

n_subs =                                length(subs);
results =                               cell( n_subs, 1 );
if classifier.do_parallel
    parfor i=1:n_subs
        results{i} =                    mrj_eeg_classify_doclassification_traintest_onesub( subs(i), feature_selection, classifier, i );
    end
else
    for i=1:n_subs
        results{i} =                    mrj_eeg_classify_doclassification_traintest_onesub( subs(i), feature_selection, classifier, i );
    end
end

%------------------------------------------------------------------------------------------------------------------
%------------------------------------------------------------------------------------------------------------------
function results = mrj_eeg_classify_doclassification_traintest_onesub( sub, feature_selection, classifier, subnum )
% unlike cross-validation, results will not be a single-item struct, but instead a 2-D array of structs
%  - results(i,j) will be struct formatted before, where i=training set and j=test set
%  - (n.b.: the above i,j are theoretical indices -- the loops below use different indices)

n_trainsets =                           length( sub.trainset );
n_testsets =                            length( sub.testset );

% initialize output struct
results(n_trainsets,n_testsets).n_features =    zeros( classifier.nits, 1 );
results(n_trainsets,n_testsets).acts =          cell( classifier.nits, 1 );
results(n_trainsets,n_testsets).testtargs =     cell( classifier.nits, 1 );

% figure out some preliminary classifier stuff
n_features =                            size( sub.trainset(1).classifier_data{1}, 1 );

for i=1:n_trainsets
    n_trainconds =                      sub.trainset(i).n_conds;
    n_traintrials_min =                 min( sub.trainset(i).n_trials_per_cond );
    traintargs =                        kron( eye(n_trainconds), ones(1,n_traintrials_min) );
    trainpats =                         zeros( n_features, n_traintrials_min * n_trainconds );
    
    for j=1:classifier.nits
        if classifier.shuffle_data_randomly
            error('Data shuffling not currently supported for separate training/test datasets');
        end
        
        for k=1:n_trainconds
            inds =                      randperm( size(sub.trainset(i).classifier_data{k},2) );
            inds =                      inds( 1:n_traintrials_min );
            trainpats(:,logical(traintargs(k,:))) = ...
                                        sub.trainset(i).classifier_data{k}(:,inds); %indexed nfeatures x ntrials
        end
        
        if feature_selection.use
            error('Feature selection not yet supported for separate training/test datasets');
        end
        % need to set results...n_features below in k loop at some point
        
        for k=1:n_testsets
            % display update
            disp([  'Subject ' int2str(subnum) ', training set ' int2str(i) ' of ' int2str(n_trainsets) ...
                    ', iteration ' int2str(j) ' of ' int2str(classifier.nits) ', test set ' int2str(k) ' of ' int2str(n_testsets)]);
            
            % train and test may not actually have the same # of conditions, so need to fake testtargs
            n_testtrials_total =        sum( sub.testset(k).n_trials_per_cond );
            n_testtrials_perfakecond =  ceil( n_testtrials_total / n_trainconds );
            testtargs_fake =            kron( eye(n_trainconds), ones(1,n_testtrials_perfakecond) );
            testtargs_fake =            testtargs_fake( :, 1:n_testtrials_total );
                                                            % # of test trials may not be evenly divided by # of fake conditions
            testpats =                  zeros( n_features, n_testtrials_total );
            
            testpat_start_ind =         1;
            for m=1:sub.testset(k).n_conds
                testpat_end_ind =       testpat_start_ind + sub.testset(k).n_trials_per_cond(m) - 1;
                testpats(:,testpat_start_ind:testpat_end_ind) = ...
                                        sub.testset(k).classifier_data{m}(:,:);
                                                            % important: testpats (and thus results) will use all test trials
                                                            %  - order will be all trials of (real) condition 1 of test set in chronological 
                                                            %    order, followed by all trials of (real) condition 2 of test set in 
                                                            %    chronological order, etc.
                testpat_start_ind =     testpat_start_ind + sub.testset(k).n_trials_per_cond(m);
            end
            
            % do actual classification here
            s =                         feval( classifier.trainfunc, trainpats, traintargs, classifier.args );
            [acts s] =                  feval( classifier.testfunc, testpats, testtargs_fake, s ); %#ok<NASGU>
            
            % don't do special behavior for SVM in this case because testtargs are fake and results will be garbage (right?)
            
            results(i,k).n_features(j) =        size(trainpats,1);
            results(i,k).acts{j} =              acts;
            results(i,k).testtargs{j} =         'testtargs invalid'; % because they're fake!
        end
    end
end
