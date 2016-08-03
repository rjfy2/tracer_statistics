function loglik = priorLoglik(mus,beta)
%prior used in sensitivity analysis: the same normal for all connections

loglik = sum(sum(log(normpdf(mus,-2,2))))+log(exppdf(beta(1),1));





