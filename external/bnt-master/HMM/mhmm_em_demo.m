% if 1
%   O = 4;
%   T = 10;
%   nex = 50;
%   M = 2;
%   Q = 3;
% else
%   O = 8;          %Number of coefficients in a vector 
%   T = 420;         %Number of vectors in a sequence 
%   nex = 1;        %Number of sequences 
%   M = 1;          %Number of mixtures 
%   Q = 6;          %Number of states 
% end
% cov_type = 'full';

% data = randn(O,T,nex);
T = 10020;   
Q = 4;
M = 1;
O = 3;
% "true" parameters
prior = normalise(rand(Q,1));
transmat = mk_stochastic(rand(Q,Q));
obsmat = mk_stochastic(rand(Q,O));
mu1 = [1 1 1]; Sigma1 = [2 .2 .1; .2 3 .2; 0.1 0.2 1];
mu2 = [2 2 2]; Sigma2 = [1 .2 .1; .2 3 .2; 0.1 0.2 4];
mu3 = [3 3 3]; Sigma3 = [3 .2 .1; .2 2 .2; 0.1 0.2 3];
mu4 = [4 4 4]; Sigma4 = [2 .2 .1; .2 2 .2; 0.1 0.2 2];
clear mu Sigma
mu(:,1) = mu1; 
mu(:,2) = mu2; 
mu(:,3) = mu3; 
mu(:,4) = mu4; 
Sigma(:,:,1) = Sigma1; 
Sigma(:,:,2) = Sigma2;
Sigma(:,:,3) = Sigma3;
Sigma(:,:,4) = Sigma4;

% training data
[obs, hidden] = mhmm_sample(T, 1, prior, transmat, mu, Sigma);

% initial guess of parameters
prior0 = normalise(rand(Q,1));
transmat0 = mk_stochastic(rand(Q,Q));

if 0
  Sigma0 = repmat(eye(O), [1 1 Q M]);
  % Initialize each mean to a random data point
  indices = randperm(T);
  mu0 = reshape(obs(:,indices(1:(Q*M))), [O Q M]);
  mixmat0 = mk_stochastic(rand(Q,M));
else
  [mu0, Sigma0] = mixgauss_init(Q*M, obs, cov_type);
  mu0 = reshape(mu0, [O Q M]);
  Sigma0 = reshape(Sigma0, [O O Q M]);
  mixmat0 = mk_stochastic(rand(Q,M));
end

[LL, prior1, transmat1, mu1, Sigma1, mixmat1] = ...
    mhmm_em(obs, prior0, transmat0, mu0, Sigma0, mixmat0, 'max_iter', 5000);


loglik = mhmm_logprob(data, prior1, transmat1, mu1, Sigma1, mixmat1);



% 
% 
% % initial guess of parameters
% prior0 = normalise(rand(Q,1));
% transmat0 = mk_stochastic(rand(Q,Q));
% 
% if 0
%   Sigma0 = repmat(eye(O), [1 1 Q M]);
%   % Initialize each mean to a random data point
%   indices = randperm(T*nex);
%   mu0 = reshape(data(:,indices(1:(Q*M))), [O Q M]);
%   mixmat0 = mk_stochastic(rand(Q,M));
% else
%   [mu0, Sigma0] = mixgauss_init(Q*M, data, cov_type);
%   mu0 = reshape(mu0, [O Q M]);
%   Sigma0 = reshape(Sigma0, [O O Q M]);
%   mixmat0 = mk_stochastic(rand(Q,M));
% end
% 
% [LL, prior1, transmat1, mu1, Sigma1, mixmat1] = ...
%     mhmm_em(data, prior0, transmat0, mu0, Sigma0, mixmat0, 'max_iter', 5);
% 
% 
% loglik = mhmm_logprob(data, prior1, transmat1, mu1, Sigma1, mixmat1);
% 
