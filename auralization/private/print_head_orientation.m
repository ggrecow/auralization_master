function print_head_orientation( spherical_angles_HRTF, head_orientation, receiver )
% function print_head_orientation( spherical_angles_HRTF, head_orientation, receiver )
%
% Prints (console and log.txt) the incidence angles of the DIRECT path used for the
% binaural rendering at the first and last time step, so that the head orientation can be
% checked in the log. For each time step, it prints
%
%   - azimuth and elevation from which the direct sound arrives at the receiver (see
%     <get_arrival_angles>; azimuth: 0 deg = +x, 90 deg = +y)
%   - HRTF azimuth = mod( arrival azimuth + head_orientation, 360 ), i.e. the head-related
%     azimuth passed to <AKhrirInterpolation> (FABIAN convention: 0 = front, 90 = left ear,
%     180 = back, 270 = right ear)
%
% INPUTS
%   spherical_angles_HRTF : struct, output of <get_propagation> (fields direct_path and
%                           reflected_path, [azimuth elevation] in degrees per time step)
%   head_orientation : scalar, in degrees (see <get_head_orientation>)
%   receiver : [x y z] receiver position, in meters
%
% Gil Felix Greco, Braunschweig 14.11.2026
%
%   AI disclosure: Generative AI assistance (Anthropic Claude) was used to develop this
%   code under the author's direction. All codes were verified by the author.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

angles = spherical_angles_HRTF.direct_path;

facing = mod( -head_orientation, 360 ); % azimuth the head faces (0 = +x, 90 = +y)

fprintf( '\n- Head orientation check (binaural rendering, direct path), receiver at x = %.5g, y = %.5g, z = %.5g m\n', ...
    receiver(1), receiver(2), receiver(3) );
fprintf( '  head_orientation = %g deg -> head faces azimuth %g deg in the x-y plane%s\n', ...
    head_orientation, facing, il_axis_label( facing ) );

idx = [1 size(angles,1)];
label = {'first', 'last'};
for k = 1:2
    az = mod( angles(idx(k),1), 360 );
    el = angles(idx(k),2);
    azHRTF = mod( az + head_orientation, 360 );
    fprintf( '  %-5s time step: sound arrives from azimuth %5.1f deg, elevation %5.1f deg -> HRTF azimuth %5.1f deg (%s)\n', ...
        label{k}, az, el, azHRTF, il_sector( azHRTF ) );
end

end

function s = il_sector( az )
% head-related sector of the HRTF azimuth (FABIAN convention)
if az < 45 || az >= 315
    s = 'front';
elseif az < 135
    s = 'left';
elseif az < 225
    s = 'back';
else
    s = 'right';
end
end

function s = il_axis_label( az )
% name of the coordinate axis for the main directions
switch az
    case 0,   s = ' (+x)';
    case 90,  s = ' (+y)';
    case 180, s = ' (-x)';
    case 270, s = ' (-y)';
    otherwise, s = '';
end
end
