%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Implementation of Kernel PCA random realization generation in MATLAB/GNU Octave
%   Author: Osho Ilamah.
%

%%  INITIALIZATION
%   Global variables/settings should be provided here
clc;
clear all;
close all;
addpath ('scripts');

%
m = 100;
n = 2500;
k = 50;
N = 50;

%
INPUT = csvread('J:\Operational\RETE\ilamaho\private\STAGE\WQ_MRW\2017\KPCA\Setup\Petrel_Examples\grid\OUT\UNPACK_MERGE.txt');
%INPUT = csvread('J:\Operational\RETE\ilamaho\private\STAGE\WQ_MRW\2017\KPCA\Setup\Petrel_Examples\grid\OUT\test.csv');
INPUT = reshape(INPUT, [n, m]);
INPUT = INPUT';

%%   Mean centering of input dataset
%
xmean =sum(INPUT, 1)/m;
XMean = repmat(xmean, m, 1);

X = INPUT - XMean;

clear INPUT XMean;

%   Specify the kernel options
%   Modify Polynomial kernel properties in this section
%      The settings here will be ignored if GRBF kernel is chosen below
kernel.order =  3; 
kernel.offset = 1; 
%
%   Modify GRBF kernel properties in this section
%      The settings here will be ignored if Polynomial kernel is chosen below
%   TODO
kernel.sigma = 7.132022; % TODO for example use the heuristics by C. Lampert to estimate sigma from X

%   Select the kernel type
kernel.type = 'polynomial';
%kernel.type = 'GRBF';

%   Specify the preimage methods
preimage_method = 'EXPLICIT';
%preimage_method = 'MDS';



%%  STANDARD PCA -Compute covariance C and decompose directly
%{
C = (1/m)*(X'*X);

% timed eigen decomposition
tic
[Qk, Dk] = eigs(C, k, 'lm');
fprintf('Direct truncated eigen decomposition took : %f\n\n', toc);

% timed svd decomposition
##tic
##[Uk, Sk, Vk] = svds(C, k);
##fprintf('Direct truncated svd decomposition took : %f\n\n', toc);
##
##Zk = X * Uk;
##
##[Qkk, Dkk] = eigs(C, 60, 'lm');
%}

%%  KERNEL PCA
%   Generate uncentered kernel matrix
Kraw = computeKernelMatrix(X, kernel);
%   Center the kernel matrix
OneM = ones(m, m)/m;
%K = Kraw - OneM*Kraw - Kraw*OneM + OneM*Kraw*OneM; %we use the more compact methos by Hai
H = eye(m) - OneM;
K = H * Kraw * H;  % equivalent to line 70 above, but more compact

%clear Kraw;
clear OneM;

%   Use eigs for full decomposition, and get ordering as bonus
tic
[Qr, Dr] = eigs(K, m, 'lm');
fprintf('Kernel truncated eigen decomposition took : %f\n\n', toc);

%{
tic
[Qr, Dr] = eig(K);
fprintf('Kernel full eigen decomposition took : %f\n\n', toc);

% Sort the eigenvalues and the eigenvectors in descending order of eigenvalues.
[Dr, decreasingOrder] = sort(Dr, 'descend'); 
Qr = Qr(:, decreasingOrder); 
%}

Dr = diag(Dr);
Dr(Dr<0) = 0;  %K is PSD and eigenvalues must all e >= 0


lambdaRm = Dr; % first m largest eigenValues of feature space centered Kernel matrix
%lambdaR = lambdaRm/m;  % first m largest eigenValues of feature space covariance matrix

sqrtLambdaRm = sqrt(lambdaRm); % sqrt of first m largest eigenValues K or first m largest (m*eigenValues) of feature space C
%   Normalize the eigenvectors of the centered Kernel matrix - we are using Dr here, c.f. PCA XXT trick and Murphy, same in essence
for ii = 1: m
    if lambdaRm(ii) ~= 0
        % normalization is only done for vectors with non-zero eigen values. Only the last eigen pair will be zero
        % but this is likely be excluded in any case
        Qr(:, ii) = Qr(:, ii) / sqrtLambdaRm(ii);
    end
end

%   At this stage:
%       Qr corresponds to A in  Hai Vo's thesis/paper (Thesis eqn 4.12], just not yet truncated
%       sqrtLambdaRm corresponds to sigma~ in  Hai Vo's thesis/paper (Thesis eqn 4.12], just not truncated

%   Now we truncate  according to the chosen k, A, SIGMA and B are defined in Hai's Thesis

A = Qr(:, 1:k);
SIGMA =  diag(sqrtLambdaRm(1:k)); % this is now a matrix square root of eigenvalues of K~ 

clear Qr Dr lambdaRm sqrtLambdaRm decreasingOrder;

%
B = H * A * (SIGMA/sqrt(m-1));


fprintf('Norm of difference in Beta matrices from 2 methods: %e \n', norm(BETA - beta_2, 'fro'));
fprintf('Norm of difference in B matrices from 2 methods: %e \n', norm(B - B_2, 'fro'));
fprintf('Norm of difference in Kraw matrices from 2 methods: %e \n', norm(Kraw_2 - Kraw, 'fro'));
fprintf('Norm of difference in Kc matrices from 2 methods: %e \n', norm(Kc_2 - K, 'fro'));


%% ***************************   *********************************   **********************************
%% Generate a new feature space realization based on a standard normal random vector ei
ei = randn(k, 1);
BETA = (B * ei) + (ones(m, 1)/m); 


clear A B;

%{ SOLVE THE PREIMAGE PROBLEM


[x, ni] = getPreimage(BETA, Kraw, X, N, preimage_method, kernel);
%[x, ni] = getPreimage(BETA, K, X, N, preimage_method, kernel);

%}

%%{
  % Optionally, plot realization 
  xi = X(ni(1), :) + xmean;
  
  x = x' + xmean;
  
  visualize2D(xi, 50, 50, 'Example realization');
  visualize2D(x, 50, 50, 'Created realization');
  % clear x xi;  
%}





