
function K = computeKernelMatrix(X, options)
% COMPUTEKERNELMATRIX  computes the kernel matrix for a centered, row organized data matrix, using the provided options   
%
%   K = COMPUTEKERNELMATRIX(X, options) computes the kernel matrix of the centered, row organized data matrix, using 
%   the chosen kernel. 
%
%   Examples: <provide examples here>
%
%   Notes:  <provide notes here>
%
%
%   See also: <list related files and functions here>
%
%   $Author: Ilamah, Osho $ $Date:2018.07.03 $ $Revision: 0.1  
 
if strcmp(options.type, 'polynomial')
    % polynomial kernel
    K = ((X*X')+options.offset).^options.order;
elseif strcmp(options.type, 'GRBF')
    % GRBF kernel
    sigma = options.sigma;
    squaredDistances = sdm(X);   
    %fprintf('SIGMA: %e \n', sigma);
    K = exp((-squaredDistances) / (2 * sigma * sigma));
    
else
    fprintf ('Unsupported kernel type');
    n = size(X, 2);
    K = zeros(n, n);
end

end