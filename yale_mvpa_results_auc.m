function auc_out = yale_mvpa_results_auc( results )
% auc_out will be n_conds x 1 matrix of AUC values (should be > .5 if better than chance)

n_conds =                               size( results.acts,1 );
auc_out =                               zeros( n_conds, 1 );

acts =                                  results.acts;
testtargs =                             results.testtargs;
    
for i=1:n_conds
    [tp, fp] =                          roc( testtargs(i,:)', acts(i,:)' );
    auc_out(i) =                        auroc( tp, fp );
end
