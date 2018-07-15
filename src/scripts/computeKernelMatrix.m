
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
    squaredDistances = squaredDistanceMatrix(X);
    
    %{
       X2 = X';
       K2 = bsxfun(@plus, dot(X2, X2, 1)', dot(X2, X2, 1)) - 2 * (X * X2);      
       fprintf('Norm of difference in Distance matrices from 2 methods: %e \n', norm(squaredDistances - K2, 'fro'));    
    %}
    
    
    fprintf('SIGMA: %e \n', sigma);
    
    K = exp((-squaredDistances) / (2 * sigma * sigma));
    
else
    fprintf ('Unsupported kernel type');
    n = size(X, 2);
    K = zeros(n, n);
end

end


function DX2 = squaredDistanceMatrix(X)
  m   = size(X, 1);
  Y   = dot(X, X, 2);
  Z1  = repmat(Y,1,m);
  Z2  = repmat(Y',m,1);
  DX2 = Z1 + Z2 - (2*(X*X'));
  % negative distances have no meaning and should not happen
  DX2(DX2 < 0) = 0 ;
end