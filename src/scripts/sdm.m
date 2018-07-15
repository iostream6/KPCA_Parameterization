
function dSquared = sdm(Xrow)
    % GETSDM  computes the "squared distances" associated with vectors of the input matrix
    % 
    %   Syntax: dSquared = sdm(Xrow)
    %
    %   Inputs: 
    %       Xrow  - The input matrix, the "squared distance" between each ROW vector constituting this matrix is to be computed
    %
    %   Outputs: 
    %   dSquared  - the "squared distances" matrix
    
    %   Examples: <provide examples here>
    %
    %   Notes:  <provide notes here>
    %
    %   References:  
    %
    %   See also: <upper case comma separated list of related files and functions here>
    %
    %   $Author: Ilamah, Osho $ $Date:2018.07.15 $ $Revision: 0.1 

      m   = size(Xrow, 1);
      Y   = dot(Xrow, Xrow, 2);
      Z1  = repmat(Y,1,m);
      Z2  = repmat(Y',m,1);
      dSquared = Z1 + Z2 - (2*(Xrow*Xrow'));
      
      % negative squared distances have no meaning and should not exist to begin with
      dSquared(dSquared < 0) = 0 ;
end
    
    