function [angles, deviation] = get_arrival_angles( ray )
% function [angles, deviation] = get_arrival_angles( ray )
%
% Direction from which the sound of an eigenray ARRIVES at the receiver, used for the HRTFs.
% It is taken from the wavefront normal n at the last point of the eigenray (at the
% receiver): the direction receiver -> sound is d = -n_end.
%
% For straight rays (homogeneous atmosphere without wind), d equals the reversed launch
% direction at the source (-n0). With refraction and wind, both differ. The launch direction
% describes the EMISSION direction (source directivity) and is not suitable for the HRTFs.
% For a ground-reflected eigenray, the last normal points upwards (after the reflection),
% so the sound arrives from below (negative elevation).
%
% INPUT
%   ray : AtmosphericRay (ART), e.g. eigenrays(i,1) (direct) or eigenrays(i,2) (reflected)
%
% OUTPUTS
%   angles : [azimuth elevation], in degrees, FABIAN convention for a head facing +x
%            azimuth in [0, 360): 0 = +x (front), 90 = +y (left), 180 = -x (back), 270 = -y (right)
%            elevation: 90 = above, 0 = horizontal, -90 = below
%   deviation : angle between d and the reversed launch direction -n0, in degrees
%               (0 for straight rays; only meaningful for rays without reflection)
%
% Gil Felix Greco, Braunschweig 14.11.2026
%
%   AI disclosure: Generative AI assistance (Anthropic Claude) was used to develop this
%   code under the author's direction. All codes were verified by the author.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nEnd = ray.n.cart(end,:);
d = -nEnd / norm( nEnd ); % unit vector, receiver -> where the sound comes from

azimuth = mod( atan2d( d(2), d(1) ), 360 );
elevation = asind( min( 1, max( -1, d(3) ) ) );
angles = [ azimuth, elevation ];

if nargout > 1
    n0 = ray.n0.cart;
    d0 = -n0 / norm( n0 ); % reversed launch direction
    deviation = acosd( min( 1, max( -1, dot( d, d0 ) ) ) );
end

end
