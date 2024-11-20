function [out] = ita_openGL2sketchUp(in)
%ITA_SU2OPENGL Summary of this function goes here
%   Detailed explanation goes here
if isa(in,'itaCoordinates')
    in=in.cart;
end
if ~(size(in,2)==3)
    error('Input has to be itaCoordinates or nx3 matrix')
end

out=itaRavenProject.pRVN2SU(in);

end

