
function [mu, sigma] = combineMeanStd(mu_i, sigma_i, N_i)

%https://en.wikipedia.org/wiki/Pooled_variance#Sample-based_statistics

assert(length(mu_i)==length(sigma_i));

N = sum(N_i);

mu = 1./N *sum( N_i.*mu_i);

term1 = sum(   (N_i-1).*sigma_i.^2 + N_i.*mu_i.^2    );
term2 = N*mu.^2;
varTotal = 1./(N-1) *(term1 -term2);
sigma = sqrt(varTotal);