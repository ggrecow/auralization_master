function doppler = get_doppler_method()
% function doppler = get_doppler_method()
%
% Reads how the Doppler effect is included in the auralization (<doppler_method> in the
% .ini file) and whether the convective amplification is applied (<convective_amplification>),
% and checks both values. The method is a decision of the user, made for each run: it
% must match the sound emission predictions (e.g. PANAM), which cannot be checked by the code.
%
%   doppler_method = input (DEFAULT)
%       The sound emission predictions already contain the Doppler effect (Doppler-shifted
%       frequencies). The signals are synthesized from these data, and the propagation
%       filters only add the propagation effects. The buzzsaw tones, which are
%       calculated from the engine speed, are Doppler-shifted in <getTonalInput> with the
%       emission angle of the input data. The time steps are placed at the receiver at constant
%       intervals (dt), i.e. the timing at the receiver follows the timing of the source; the
%       Doppler shift of the frequencies comes from the predictions.
%       <convective_amplification> is not used: if needed, it must be in the input data.
%
%   doppler_method = delay_line
%       The sound emission predictions do NOT contain the Doppler effect. The source signal is
%       synthesized in emission time. Each time step is placed at the time it reaches the
%       receiver (emission time + propagation time of the eigenray), separately for the
%       direct and ground-reflected sound paths, i.e. at uneven intervals (closer together
%       while the aircraft approaches, further apart while it moves away). Compressing and
%       stretching the signal between the steps produces the Doppler shift (variable delay
%       lines, see <get_doppler_mapping> and <apply_doppler_delay_line>), so pitch and timing
%       at the receiver both follow from the propagation times. Using Doppler-shifted input
%       data with this method applies the Doppler effect twice.
%
%   convective_amplification = 0 (DEFAULT) or 1, only used with doppler_method = delay_line
%       0 -> no convective amplification (e.g. input data already contain it).
%       1 -> the level of each sound path is changed by -40*log10(D), with the Doppler
%            factor D = 1 - M*cos(theta) (+5.0 dB for M = 0.25 straight towards the receiver).
%            Use this if the sound emission predictions do not contain it.
%
% OUTPUT : struct
%   doppler.method                   : 'input' or 'delay_line'
%   doppler.convective_amplification : 0 or 1
%
% Gil Felix Greco, Braunschweig 16.11.2026
%
%   AI disclosure: Generative AI assistance (Anthropic Claude) was used to develop this
%   code under the author's direction. All codes were verified by the author.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global input_file

% doppler_method
if isfield( input_file, 'doppler_method' )
    method = lower( strtrim( char( input_file.doppler_method ) ) );
    if ~any( strcmp( method, {'input', 'delay_line'} ) )
        error( 'auralization:invalidDopplerMethod', ...
            '<doppler_method> in the .ini file must be input or delay_line, found: %s', method );
    end
else
    method = 'input'; % default value
end

% convective_amplification
if isfield( input_file, 'convective_amplification' )
    convective_amplification = str2double( input_file.convective_amplification );
    if ~( isequal( convective_amplification, 0 ) || isequal( convective_amplification, 1 ) )
        error( 'auralization:invalidConvectiveAmplification', ...
            '<convective_amplification> in the .ini file must be 0 or 1, found: %s', ...
            strtrim( char( input_file.convective_amplification ) ) );
    end
else
    convective_amplification = 0; % default value
end

doppler.method = method;
doppler.convective_amplification = convective_amplification;

end
