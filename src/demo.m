% demo for running BCD-IC
% add the MEX files to the path
addpath('MEX');
addpath('clustering\metis-5.0.2\metismex');

%%%%%%%%%%%%%%%%%%%% GENERATE THE PROBLEM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% problem dimensions
n = 2048; % dimension

k = 100; % number of samples
lambda = 0.42;
% Create the matrix X
rng(1); % set random seed for creating always the same example
[A_true,X] = generateSyntheticExperiment(n,k,'random_planar_Lap',0.01);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Set the parameters for BCD-IC
params = [];
params.epsilon = 1e-2; % this is epsilon for the stopping criterion 
params.epsilon_threshold = 1e-12;
params.verbose = true;
params.max_iter = 20;
params.CGiter = 100; % maximum number of iterations for LASSO (NLCG)
% params.blockSize = 512;
params.blockSize = 10;

params.epsilonMGi = 1e-4;
params.epsilonInvAi = 1e-4;
params.epsilonInnerCG = 1e-4;
params.num_cores = 4;
params.LinearSolution = 'CG';
params.MultilevelAcceleration = false;
params.devideAndConquer = false;
params.continuation = false;
params.PCG_blockSize = 2^14;

% Run BCD-IC

X = X'; % we use the data as a mXn matrix. 
params.MultilevelAcceleration = true;
[A,f,support,active,timeSamples] = BCDIC(X,lambda,params);
params.MultilevelAcceleration = false;
[A,f,support,active,timeSamples] = BCDIC(X,lambda,params);
