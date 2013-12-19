function [yale_mvpa_config, subs, results] = yale_mvpa_doclassification( yale_mvpa_config, subs )

if ~ischar(yale_mvpa_config.classifier.kfold) || strcmp(yale_mvpa_config.classifier.kfold,'subs')
    error('not supported yet');
end

if ~isempty(yale_mvpa_config.classifier.nits) && yale_mvpa_config.classifier.nits ~= 0
    error('not supported yet');
end

if yale_mvpa_config.classifier.shuffle_data_randomly ~= 0
    error('not supported yet');
end

% the above should morph into checks for mutually incompatible options at some point when features are actually put back in

yale_mvpa_config.results.function =     @yale_mvpa_display_process_results; % at some point, presumably we'll want to choose something else...
                                                                            %  likely due to combinations of the options above

n_subs =                                length(subs);
results =                               cell( n_subs, 1 );
if yale_mvpa_config.general.try_parallel
    parfor i=1:n_subs
        results{i} =                    yale_mvpa_doclassification_onesub( subs, yale_mvpa_config, i );
    end
else
    for i=1:n_subs
        results{i} =                    yale_mvpa_doclassification_onesub( subs, yale_mvpa_config, i );
    end
end




%-------------------------------------------------------------------------------------
function results = yale_mvpa_doclassification_onesub( subs, yale_mvpa_config, subnum )

nits =                                                      1; %at some point, change to reflect yale_mvpa_config.classifier.nits
results =                                                   struct('n_features',cell(1,nits),'acts',cell(1,nits),'testtargs',cell(1,nits));

for i = 1:nits
    % display update
    fprintf('Subject %d, iteration %d of %d\n', subnum, i, nits);
    
    [trainpats, traintargs, testpats, testtargs] =          yale_mvpa_doclassification_divide_traintest( subs, yale_mvpa_config, subnum );
    
    if yale_mvpa_config.classifier.mean_treatment == 1 %subtract out mean
        for j = 1:numel(trainpats)
            trainpats{j} = trainpats{j} - repmat( mean(trainpats{j}), size(trainpats{j},1), 1);
            testpats{j}  = testpats{j}  - repmat( mean(testpats{j}),  size(testpats{j},1),  1);
        end
    elseif yale_mvpa_config.classifier.mean_treatment == 2 %keep only mean
        for j = 1:numel(trainpats)
            trainpats{j} = mean(trainpats{j});
            testpats{j}  = mean(testpats{j});
        end
    elseif yale_mvpa_config.classifier.mean_treatment == 3 %z-score
        for j = 1:numel(trainpats)
            trainpats{j} = zscore(trainpats{j});
            testpats{j}  = zscore(testpats{j});
        end
    end
    
    %some form of feature selection would probably happen here
    
    n_folds =                                               length(trainpats); %folds, runs, what-have-you
    n_conds =                                               size(traintargs{1},1); %could probably move out of loop, but meh
    n_test_trials_per_fold =                                size(testtargs{1},2); %assume for now it is the same for all runs/folds
    acts =                                                  nan(n_conds,n_test_trials_per_fold,n_folds);
    all_testtargs =                                         acts;
    
    for k = 1:n_folds
        class_struct =                                      feval( yale_mvpa_config.classifier.trainfunc, trainpats{k}, traintargs{k}, yale_mvpa_config.classifier.args );
        [acts(:,:,k), class_struct] =                       feval( yale_mvpa_config.classifier.testfunc,  testpats{k},  testtargs{k},  class_struct ); %#ok<NASGU>
        all_testtargs(:,:,k) =                              testtargs{k};
    end
    
    results(i).n_features =                                 size(trainpats{1},1);
    results(i).acts =                                       reshape(acts,[n_conds numel(acts)/n_conds]);
    results(i).testtargs =                                  reshape(all_testtargs,[n_conds numel(all_testtargs)/n_conds]);
end