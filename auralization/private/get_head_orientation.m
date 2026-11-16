function head_orientation = get_head_orientation()
% function head_orientation = get_head_orientation()
%
% Reads the head orientation used for the binaural rendering (<head_orientation> in the
% .ini file, in degrees) and checks that it is a finite number. DEFAULT value: 0.
%
% <head_orientation> is a fixed rotation (yaw) of the head in the horizontal (x-y) plane of
% the PANAM coordinate system. It is the same for all receivers and does not follow the
% aircraft. The HRIRs are interpolated at the azimuth
%
%   (azimuth of the source seen from the receiver) + head_orientation
%
% where the source azimuth is 0 deg = +x and 90 deg = +y. Increasing values therefore turn
% the head CLOCKWISE seen from above (head facing direction = -head_orientation):
%   head_orientation = 0   -> head faces +x (DEFAULT)
%   head_orientation = 90  -> head faces -y
%   head_orientation = 180 -> head faces -x
%   head_orientation = 270 -> head faces +y
% See the header of <get_FIR> for what the listener hears in each case.
%
% OUTPUT
%   head_orientation : scalar, in degrees
%
% Gil Felix Greco, Braunschweig 14.11.2026
%
%   AI disclosure: Generative AI assistance (Anthropic Claude) was used to develop this
%   code under the author's direction. All codes were verified by the author.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global input_file

if isfield( input_file, 'head_orientation' )
    head_orientation = str2double( input_file.head_orientation ); % get <head_orientation> from <input_file>
    if ~isscalar( head_orientation ) || ~isfinite( head_orientation )
        error( 'auralization:invalidHeadOrientation', ...
            '<head_orientation> in the .ini file must be a number (in degrees), found: %s', ...
            strtrim( char( input_file.head_orientation ) ) );
    end
else
    head_orientation = 0; % default value
end

end
