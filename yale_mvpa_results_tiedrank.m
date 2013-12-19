function tiedrank_out = yale_mvpa_results_tiedrank( results )
% tiedrank should be 1 x ntrials, containing rank score for the correct category on each trial
%  (scale is 1-n_conds, with higher being better; chance is n_conds/2)

acts =                              results.acts;
testtargs =                         logical(results.testtargs);
tiedrank_tmp =                      tiedrank(acts);
tiedrank_out =                      tiedrank_tmp(testtargs);
