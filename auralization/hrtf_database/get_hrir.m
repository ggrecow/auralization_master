function [hrirLeft, hrirRight] = get_hrir( azimuth, elevation, hrtf )
% function [hrirLeft, hrirRight] = get_hrir( azimuth, elevation, hrtf )
%
% Head-related impulse responses (HRIRs) of the HRTF set <hrtf> (see <get_hrtf_database>)
% for the given directions, interpolated from the spherical-harmonics (SH) representation
% of the MEASURED HRTFs (HATO = 0):
%
%   FABIAN : AKhrirInterpolation( azimuth, elevation, 0, 'measured_sh' )
%   HUTUBS : AKhrirDatabase( <SH coefficients>, [azimuth elevation], false, 'mat' ), i.e.
%            the 'sh' mode of AKhrirDatabase for measured data (SH interpolation of the
%            time-aligned HRTFs, then re-insertion of the interpolated times of arrival).
%            The coefficients are passed as a struct, so the file is not read again.
%
% Coordinates (head-centred, the same for both databases): azimuth 0 = front, 90 = left
% ear, 180 = back, 270 = right ear; elevation 90 = above, 0 = front, -90 = below. Azimuths
% outside [0, 360) are allowed.
%
% INPUTS
%   azimuth, elevation : [Q x 1] (or [1 x Q]) directions, in degrees
%   hrtf               : struct from <get_hrtf_database> (DEFAULT: get_hrtf_database())
%
% OUTPUTS
%   hrirLeft, hrirRight : [256 x Q] HRIRs at hrtf.fs (44.1 kHz)
%
% Gil Felix Greco, Braunschweig 19.11.2026
%
%   AI disclosure: Generative AI assistance (Anthropic Claude) was used to develop this
%   code under the author's direction. All codes were verified by the author.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3 || isempty( hrtf )
    hrtf = get_hrtf_database();
end

azimuth   = azimuth(:);
elevation = elevation(:);

if numel( azimuth ) ~= numel( elevation )
    error( 'get_hrir:input', 'azimuth and elevation must have the same number of elements.' );
end

switch hrtf.database

    case 'FABIAN'

        HATO = 0; % head-above-torso orientation (only HATO = 0 is used in the framework)
        [hrirLeft, hrirRight] = AKhrirInterpolation( azimuth, elevation, HATO, 'measured_sh' );

    case 'HUTUBS'

        hrir = AKhrirDatabase( hrtf.sh, [azimuth elevation], false, 'mat' ); % [256 x Q x 2]
        hrirLeft  = hrir(:,:,1);
        hrirRight = hrir(:,:,2);

    otherwise
        error( 'get_hrir:database', 'Unknown HRTF database ''%s''.', hrtf.database );
end

end
