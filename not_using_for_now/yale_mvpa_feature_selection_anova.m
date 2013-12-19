function inds = mrj_feature_selection_anova( trainpats, traintargs, args )
% does simple ANOVA over training data
%  takes in pre-selected training data, so assumes # of trials is the same across all conditions
% trainpats is nfeatures x ntrials
% traintargs is nconds x ntrials
% args: one field, args.thresh, indicating p-value cutoff
% returns: inds, indices of the features that pass the cutoff

% pre-allocate p-value vector, set up grouping variable
anova_ps =                              zeros(size(trainpats,1),1); % one p-val for each feature
anova_grouping =                        (traintargs') * ((1:size(traintargs,1))');

% loop over features, do ANOVA on each one
if args.parallel
    parfor i=1:size(trainpats,1)
        anova_ps(i) =                       anova1( trainpats(i,:), anova_grouping, 'off' );
    end
else
    for i=1:size(trainpats,1)
        anova_ps(i) =                       anova1( trainpats(i,:), anova_grouping, 'off' );
    end
end

inds =                                  find(anova_ps < args.thresh);
