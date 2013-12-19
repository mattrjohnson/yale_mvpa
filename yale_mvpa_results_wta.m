function wta = yale_mvpa_results_wta( results )
% wta will be 1 x ntrials matrix of 1's (winner == correct) and 0's
%  (chance is 1/n_conds, if you average all the 1's and 0's together)

acts =                              results.acts;
testtargs =                         results.testtargs;
acts_maxes =                        max(acts);
acts_tmp =                          zeros(size(acts));
for j=1:length(acts_maxes)
    acts_tmp(:,j) =                 (acts(:,j)==acts_maxes(j));
end
acts_tmp2 =                         (acts_tmp & testtargs);
wta =                               sum(acts_tmp2);
wta(wta>0) =                        1./sum(acts_tmp(:,wta>0)); %adjusts for ties -- trials with ties get fractional accuracies
