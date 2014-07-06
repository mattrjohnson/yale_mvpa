function [yale_mvpa_config, subs, results] = yale_mvpa_doclassification( yale_mvpa_config, subs )

if ischar(yale_mvpa_config.classifier.kfold) && (~isempty(yale_mvpa_config.classifier.kfold)) && (~strcmp(yale_mvpa_config.classifier.kfold,'runs'))
    error('not supported yet');
end

if isnumeric(yale_mvpa_config.classifier.kfold)
    if mod(yale_mvpa_config.classifier.kfold,1) || ~(yale_mvpa_config.classifier.kfold > 1)
        error('not supported yet, or weird numeric value for yale_mvpa_config.classifier.kfold');
    end
end

%there could be other bad combinations of cases for kfold... maybe check for them here at some later point, although
% hopefully they should be caught by the catch-all error case in yale_mvpa_doclassification_divide_traintest

% if ~isempty(yale_mvpa_config.classifier.nits) && yale_mvpa_config.classifier.nits ~= 0
%     error('not supported yet');
% end
% not sure if we need to do much else to enable multiple-iterations support, but I guess we'll find out...

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

nits =                                                      yale_mvpa_config.classifier.nits;
if isempty(nits) || (isnumeric(nits) && isscalar(nits) && nits==0)
    nits = 1;
elseif ~isnumeric(nits) || ~isscalar(nits) || nits<1 || mod(nits,1)
    error('Seemingly invalid value of yale_mvpa_config.classifier.nits!');
end

results =                                                   struct('n_features',cell(1,nits),'acts',cell(1,nits),'testtargs',cell(1,nits));

fprintf('Subject %.3d: ', subnum);
for i = 1:nits
    % display update
    update1_str = sprintf('  iteration:%.4d  of:%.4d', i, nits);
    fprintf(update1_str);
    
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
    
%     if yale_mvpa_config.classifier.cheat %add in shuffle data randomly options around here later
%         
%     end
    
    n_folds =                                               length(trainpats); %folds, runs, what-have-you
    n_conds =                                               size(traintargs{1},1); %could probably move out of loop, but meh
    n_test_trials_per_fold =                                size(testtargs{1},2); %assume for now it is the same for all runs/folds
    acts =                                                  nan(n_conds,n_test_trials_per_fold,n_folds);
    all_testtargs =                                         acts;
    
    for k = 1:n_folds
        update2_str = sprintf('  fold:%.3d', k);
        fprintf(update2_str);
        class_struct =                                      feval( yale_mvpa_config.classifier.trainfunc, trainpats{k}, traintargs{k}, yale_mvpa_config.classifier.args );
        [acts(:,:,k), class_struct] =                       feval( yale_mvpa_config.classifier.testfunc,  testpats{k},  testtargs{k},  class_struct ); %#ok<NASGU>
        all_testtargs(:,:,k) =                              testtargs{k};
        if (i ~= nits) || (k ~= n_folds)
            fprintf(repmat('\b',1,length(update2_str)));
        end
    end
    
    if i == nits
        fprintf('\n');
    else
        fprintf(repmat('\b',1,length(update1_str)));
    end
    
    results(i).n_features =                                 size(trainpats{1},1);
    results(i).acts =                                       reshape(acts,[n_conds numel(acts)/n_conds]);
    results(i).testtargs =                                  reshape(all_testtargs,[n_conds numel(all_testtargs)/n_conds]);
end