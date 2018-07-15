function visualize2D(z, NX, NY, str_title)
% VISUALIZE2D  plots a vector in 2 dimensions   
% 
%   Syntax: visualize2D(z, NX, NY, str_title)
%
%   Inputs: 
%        z    - The vector to be visualized
%        NX   - The number of cells/elements in the X direction
%        NY   - The number of cells/elements in the Y direction
%   str_title - The title of the 2D visualization
%
%
%   Notes:  <provide notes here>
%
%   See also: <upper case comma separated list of related files and functions here>
%
%   $Author: Ilamah, Osho $ $Date:2018.07.13 $ $Revision: 0.1  

%

range = [min(z), max(z)];

figure
grid_data = reshape(z, NX, NY);
% reshape traverses columns first, so we transpose to give row first travesal
grid_data = grid_data';

% display image with scaled colors
imagesc(grid_data);
colorbar;
caxis(range);
set(gca, 'ydir', 'normal');
axis equal tight
%
title(str_title);

end