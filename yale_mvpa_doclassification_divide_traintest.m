function [trainpats, traintargs, testpats, testtargs] = yale_mvpa_doclassification_divide_traintest( subs, yale_mvpa_config, subnum )

this_sub =                                                  subs(subnum);

if ischar(yale_mvpa_config.classifier.kfold) && strcmp(yale_mvpa_config.classifier.kfold,'runs')
    % remember to move relevant lines out of this 'if' statement as we fill in the other options
    
    % set up some initial values
    n_conds =                                               size(this_sub.classifier_data,1);
    n_runs =                                                size(this_sub.classifier_data,2);
    n_features =                                            size(this_sub.classifier_data{1},1);
    n_trials_percond_perrun =                               size(this_sub.classifier_data{1},2); %at some point, add a check to confirm that all runs/conditions have same numbers of trials
    n_trials_perrun =                                       n_trials_percond_perrun * n_conds;
    n_testtrials_percond =                                  n_trials_percond_perrun; %redundant in 'runs' case, but matches what we'll do in other cases
    
    % figure out traintargs and testtargs
    %  n.b.: in this case ('runs'), both traintargs and testtargs are the same across all runs, but that may not always be the case
    testtargs_tmp =                                         kron( eye(n_conds), ones(1,n_testtrials_percond) );
    traintargs_tmp =                                        nan( n_conds,n_trials_perrun,(n_runs-1) );
    for i = 1:(n_runs-1)
        traintargs_tmp(:,:,i) =                             kron( eye(n_conds), ones(1,n_trials_percond_perrun) );
    end
    traintargs_tmp =                                        reshape( traintargs_tmp, [ n_conds  n_trials_perrun*(n_runs-1) ] );
    traintargs =                                            cell(n_runs,1);
    testtargs =                                             cell(n_runs,1);
    for i = 1:n_runs
        traintargs{i} =                                     traintargs_tmp;
        testtargs{i} =                                      testtargs_tmp;
    end
    
    % finally, figure out trainpats and testpats
    
    % first, reshape the data a bit so it's easier to deal with
    classifier_data_reshaped =          nan(n_conds, n_runs, n_features, n_trials_percond_perrun);
    for i = 1:n_conds
        for j = 1:n_runs
            classifier_data_reshaped(i,j,:,:) =             this_sub.classifier_data{i,j};
        end
    end
    classifier_data_reshaped =                              permute(classifier_data_reshaped,[2 3 4 1]); %now runs x features x trial_within_cond x cond
    classifier_data_reshaped =                              reshape(classifier_data_reshaped,[n_runs n_features n_conds*n_trials_percond_perrun]); %now runs x features x trial_within_run
    
    % now, loop through runs and make our trainpats/testpats
    trainpats =                                             cell(n_runs,1);
    testpats =                                              cell(n_runs,1);
    for j = 1:n_runs
        train_runs =                                        setdiff( 1:n_runs, j );
        trainpats_tmp =                                     nan( n_features, n_conds*n_trials_percond_perrun, n_runs-1 ); %indexed features x trials_perrun x run
        for k=1:length(train_runs)
            trainpats_tmp(:,:,k) =                          squeeze( classifier_data_reshaped(train_runs(k),:,:) );
        end
        trainpats_tmp =                                     reshape( trainpats_tmp, [n_features  n_conds*n_trials_percond_perrun*(n_runs-1) ] );
        testpats_tmp =                                      squeeze( classifier_data_reshaped( j,:,: ) );
        
        trainpats{j} =                                      trainpats_tmp;
        testpats{j} =                                       testpats_tmp;
    end
elseif yale_mvpa_config.classifier.separate_traintest == 1
    %going to take a lot of shortcuts for Img_Perc now... come back later to generalize
    % - notably, for now we'll only support a single trainset / testset
    % - also, have to assume equal numbers of trials across conditions, and nits == 1
    
    % set up some initial values
    n_features =                                            size( this_sub.trainset.classifier_data{1}, 1 );
    n_trainconds =                                          this_sub.trainset.n_conds;
    n_traintrials_percond =                                 size(this_sub.trainset.classifier_data{1},2);
    n_testconds =                                           this_sub.testset.n_conds; %should really be the same # as trainconds, but worry about that later
    n_testtrials_percond =                                  size(this_sub.testset.classifier_data{1},2);
    
    traintargs_tmp =                                        kron( eye(n_trainconds), ones(1,n_traintrials_percond) );
    testtargs_tmp =                                         kron( eye(n_testconds),  ones(1,n_testtrials_percond) );
    
    trainpats_tmp =                                         zeros( n_features, n_traintrials_percond * n_trainconds );
    testpats_tmp =                                          zeros( n_features, n_testtrials_percond *  n_testconds );
    
    for i = 1:n_trainconds
        trainpats_tmp(:,logical(traintargs_tmp(i,:))) =     this_sub.trainset.classifier_data{i};
    end
    
    for i = 1:n_testconds
        testpats_tmp(:,logical(testtargs_tmp(i,:))) =       this_sub.testset.classifier_data{i};
    end
    
    traintargs =                                            { traintargs_tmp };
    testtargs =                                             { testtargs_tmp };
    trainpats =                                             { trainpats_tmp };
    testpats =                                              { testpats_tmp };
else
    error('Not supported yet');
end










%traintargs will be arranged a bit differently for 'runs' case; will contain sections by run instead of having each condition contiguous
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
