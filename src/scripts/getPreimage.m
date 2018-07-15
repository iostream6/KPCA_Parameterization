
function [x, ni] = getPreimage(BETA, K, Xrow, N, method, kernel_options)
% GETPREIMAGE  solves the kernel PCA preimage problem using either the MDS distance-based or the explicit expression methods   
% The preimage methods here are based on the description provided in the Standford University Phd Thesis by Hai Xuan Vo:
% "New Geological Parameterization for History Matching Complex Models"
% 
%   Syntax: [x, ni] = getPreimage(BETA, K, Xrow, N, method, kernel_options)
%
%   Inputs: 
%      BETA   - defined by equation 4.14 in thesis by Hai Vo
%        K    - The centred kernel matrix
%       Xrow  - The complete training samples/realizations, row based
%        N    - The number of nearest neighbors
%     method  - character string specifying the solution method. ['MDS'|'EXPLICIT']
%    options  - struct specifying the kernel  options
%
%   Outputs: 
%        x    - the calculated preimage
%        ni   - vector whose elements are the indices of the N nearest neighbors of x

%   Examples: <provide examples here>
%
%   Notes:  <provide notes here>
%
%   References:  
%      "New Geological Parameterizations for History Matching Complex Models". Hai Vo. Stanford University PhD thesis. 2015.
%      "The pre-image problem in kernel methods". James Kwok & Ivor Tsang. ICML-2003 conference paper. 2003.
%      "Multiobjective history matching with dominance and decomposition". Osho Ilamah. Journal paper (in press).
%
%   See also: <upper case comma separated list of related files and functions here>
%
%   $Author: Ilamah, Osho $ $Date:2018.07.13 $ $Revision: 0.1  

m = size(Xrow, 1);

Xcol = Xrow'; % we need X organized column wise, so transpose a copy of the input matrix

KBETA = K* BETA;
BETA_T_KBETA = BETA' * KBETA;

if strcmp(method, 'MDS')
   if strcmp(kernel_options.type, 'polynomial')
      % polynomial kernel   
      diSquared = getPolynomialKernelSquaredInputDistance(Xrow, kernel_options, KBETA, BETA_T_KBETA, m);
   elseif strcmp(kernel_options.type, 'GRBF')
      % GRBF kernel
      % compute the feature space squared distances 
      % dfSquared = 1 + ((phi_xhat)'*(phi_xhat)) - (2*(phi_xhat)'*(phi_xi))
      dfSquared = repmat(1 + BETA_T_KBETA, m, 1) - (2 * KBETA);
      diSquared = -2 * (kernel_options.sigma^2) * log(1 - (0.5*dfSquared));
      
      clear dfSquared;      
   else
      ni = 0;
      x = 0;
      return;
   end  
   [~, distanceOrder] = sort(diSquared);  
   % select data for N nearest neighbors
   diSquared_N = diSquared(distanceOrder(1:N));
   X_N = Xcol(:, distanceOrder(1:N));
   H_N = eye(N, N) - (ones(N,N)/N);
   XN_HN = X_N * H_N;
   
   [U_N, S_N, V_N] = svd(XN_HN);
   Z_N = S_N * V_N';
   
   doSquared_N =  sum(Z_N .^ 2)';
   zhat = -0.5 * pinv(Z_N') * (diSquared_N - doSquared_N);
   
   x = U_N * zhat + mean(X_N, 2);
   ni = distanceOrder(1:N);
   
   clear zhat;
   clear diSquared diSquared_N doSquared_N;
   clear X_N H_N XN_HN Z_N;
   clear U_N S_N V_N;
  
else 
   if strcmp(kernel_options.type, 'polynomial')
      % polynomial kernel     
      diSquared = getPolynomialKernelSquaredInputDistance(Xrow, kernel_options, KBETA, BETA_T_KBETA, m);
      [~, distanceOrder] = sort(diSquared);  
      % select data for N nearest neighbors
      BETA_N = BETA(distanceOrder(1:N));
      KBETA_N = KBETA(distanceOrder(1:N));
      X_N = Xcol(:, distanceOrder(1:N));
      %
      theta_numerator = BETA_N .* nthroot((KBETA_N / BETA_T_KBETA).^(kernel_options.order -1),kernel_options.order);
      
      clear KBETA_N;
      clear diSquared;
      
   elseif strcmp(kernel_options.type, 'GRBF')
      % GRBF kernel
      % compute the feature space squared distances
      % dfSquared = 1 + ((phi_xhat)'*(phi_xhat)) - (2*(phi_xhat)'*(phi_xi))
      dfSquared = repmat(1 + BETA_T_KBETA, m, 1) - (2 * KBETA);
      [~, distanceOrder] = sort(dfSquared);  
      % select data for N nearest neighbors
      dfSquared_N = dfSquared(distanceOrder(1:N));
      X_N = Xcol(:, distanceOrder(1:N));
      BETA_N = BETA(distanceOrder(1:N));
      %
      theta_numerator = BETA_N .* (1 - (0.5* dfSquared_N));

      clear dfSquared;
      clear dfSquared_N;
      
   else
      fprintf ('Unsupported kernel type');
      ni = 0;
      x = 0;
      return;
   end
   
   clear BETA_N;
  
   theta = theta_numerator/ (sum(theta_numerator));
   x = X_N * theta;
   ni = distanceOrder(1:N);
   
   
   clear theta;
   clear theta_numerator;
   clear X_N;
  
end

   clear BETA_T_KBETA;
   clear KBETA;
   clear Xcol;

end

function diSquared = getPolynomialKernelSquaredInputDistance(Xrow, kernel_options, KBETA, BETA_T_KBETA, m)
      % make sure order is odd
      if mod(kernel_options.order, 2) == 0
        kernel_options.order = kernel_options.order + 1;
      end
      
      BETA_T_KBETA_MAT = repmat(BETA_T_KBETA, m, 1);
      %
      % compute (xi)'*(xi) == xi dot xi, i.e. eqn 13c in notes
      xi_dot_xi = dot(Xrow, Xrow, 2);
      % compute (xhat)'*(xhat) == xhat dot xhat, i.e. eqn 13c in  notes
      xhat_dot_xhat = nthroot(BETA_T_KBETA_MAT, kernel_options.order) - kernel_options.offset;
      % compute (xhat)'*(xi) == xhat dot xi, i.e. eqn 13c in notes
      xhat_dot_xi = nthroot(KBETA, kernel_options.order) - kernel_options.offset;
      %
      % compute the input space squared distances 
      diSquared = xi_dot_xi + xhat_dot_xhat - (2*xhat_dot_xi);
      
      clear BETA_T_KBETA_MAT;
      clear xi_dot_xi;
      clear xhat_dot_xhat;
      clear xhat_dot_xi;
end