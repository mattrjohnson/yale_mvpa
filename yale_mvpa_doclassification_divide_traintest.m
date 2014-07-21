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
    %  n.b.: in this case ('runs'), both traintargs and testtargs are the same across all runs, but that may not always be the case (?)
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
elseif isnumeric(yale_mvpa_config.classifier.kfold) && ~mod(yale_mvpa_config.classifier.kfold,1) && (yale_mvpa_config.classifier.kfold>1)
    % divide up into K folds
    n_folds =                                               yale_mvpa_config.classifier.kfold;
    n_conds =                                               size(this_sub.classifier_data,1);
    n_features =                                            size(this_sub.classifier_data{1},1);
    n_trials_percond =                                      zeros(n_conds,1);
    for i = 1:n_conds
        n_trials_percond(i) =                               size(this_sub.classifier_data{i},2);
    end
    min_n_trials_percond =                                  min(n_trials_percond);
    n_testtrials_percond =                                  floor(min_n_trials_percond ./ n_folds);
    if n_testtrials_percond<1
        error('Not enough trials for k-fold cross validation.');
    end
    n_traintrials_percond =                                 n_testtrials_percond * (n_folds-1);
    
    % figure out traintargs and testtargs
    %  n.b.: in this case, both traintargs and testtargs are the same across all folds, but that may not always be the case (?)
    testtargs_tmp =                                         kron( eye(n_conds), ones(1,n_testtrials_percond) );
    traintargs_tmp =                                        kron( eye(n_conds), ones(1,n_traintrials_percond) );
    traintargs =                                            cell(n_folds,1);
    testtargs =                                             cell(n_folds,1);
    for i = 1:n_folds
        traintargs{i} =                                     traintargs_tmp;
        testtargs{i} =                                      testtargs_tmp;
    end
    
    % first, reshape the data a bit so it's easier to deal with
    classifier_data_reshaped =                              cell(n_conds,1);
    for i = 1:n_conds
        classifier_data_temp =                              this_sub.classifier_data{i};
        if yale_mvpa_config.classifier.cheat %very temporary! only for internal testing
%             classifier_data_temp = ones(size(classifier_data_temp)) .* i + rand(size(classifier_data_temp)) .* yale_mvpa_config.classifier.cheat_noise_factor;
            classifier_data_temp = ones(600,size(classifier_data_temp,2)) .* i + rand(600,size(classifier_data_temp,2)) .* yale_mvpa_config.classifier.cheat_noise_factor;
            n_features = 600;
        end
        trial_inds =                                        randperm(n_trials_percond(i));
        classifier_data_temp =                              classifier_data_temp(:,trial_inds(1:(n_testtrials_percond*n_folds)));
        classifier_data_reshaped{i} =                       classifier_data_temp; %still features x trials, but randomized & truncated to right size
    end
    
    % now, loop through runs and make our trainpats/testpats
    trainpats =                                             cell(n_folds,1);
    testpats =                                              cell(n_folds,1);
    for i = 1:n_folds
        train_folds =                                       setdiff(1:n_folds,i);
        trainpats_tmp =                                     nan(n_features, n_traintrials_percond, n_conds);
        testpats_tmp =                                      nan(n_features, n_testtrials_percond, n_conds );
        
        for j = 1:n_conds
            classifier_data_temp =                          classifier_data_reshaped{j};
            classifier_data_temp =                          reshape(classifier_data_temp, [n_features n_folds n_testtrials_percond]); %now features x folds x trials
            testpats_tmp(:,:,j) =                           squeeze(classifier_data_temp(:,i,:));
            classifier_data_temp =                          classifier_data_temp(:,train_folds,:);
            trainpats_tmp(:,:,j) =                          reshape(classifier_data_temp, [n_features n_traintrials_percond]);
        end
        
        trainpats{i} =                                      reshape(trainpats_tmp, [n_features n_traintrials_percond * n_conds]);
        testpats{i}  =                                      reshape(testpats_tmp,  [n_features n_testtrials_percond  * n_conds]);
    end
elseif ischar(yale_mvpa_config.classifier.kfold) && strcmp(yale_mvpa_config.classifier.kfold,'custom')
    [trainpats, traintargs, testpats, testtargs] = feval(yale_mvpa_config.classifier.kfold_custom, this_sub, yale_mvpa_config);
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
