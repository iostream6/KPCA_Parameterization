
function sigma = getSigma(Xrow, param)
   % GETSIGMA  estimates sigma for the GRBF kernel based on heuristics   
   % This code is based on MATLAB implementations provided by Hai Vo
   % 
   %   Syntax: sigma = getSigma(Xrow, param)
   %
   %   Inputs: 
   %       Xrow  - The complete training samples/realizations, row based
   %      param  - user supplied input parameter used in the estimation of sigma
   %
   %   Outputs: 
   %        sigma    - the estimated sigma for the input dataset
   
   %   Examples: <provide examples here>
   %
   %   Notes:  <provide notes here>
   %
   %   References:  
   %      "New Geological Parameterizations for History Matching Complex Models". Hai Vo. Stanford University PhD thesis. 2015.
   %      "Kernel Methods in Computer Vision". Christoph H. Lampert. http://mr.crossref.org/iPage?doi=10.1561%2F0600000027
   %      "Multiobjective history matching with dominance and decomposition". Osho Ilamah. Journal paper (in press).
   %
   %   See also: <upper case comma separated list of related files and functions here>
   %
   %   $Author: Ilamah, Osho $ $Date:2018.07.15 $ $Revision: 0.1  

   squaredDistances = sdm(Xrow);
   % unrolling into a vector
   squaredDistances = squaredDistances(:);
   % Modified from Lampert, following Vo
   sigma = param * sqrt(median(squaredDistances) / 2);

end