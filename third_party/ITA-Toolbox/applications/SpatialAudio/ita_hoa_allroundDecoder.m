function [D_allrad] = ita_hoa_allroundDecoder(LS_pos, tDesign,varargin)
%ALLRAD Calculate an Allround Ambisonics Decoder (AllRAD)
% The AllRAD uses a concatenation of Ambisonics and VBAP to overcome
% reproduction setup irregularities such as inconsistent spatial coverage
% and non-spherical setups. This is done by means of a virtual loudspeaker
% array ideally suited for HOA reproduction whose output channels are
% VBAPed onto the real loudspeaker channels.
% Caution: Converges toward VBAP high Ambisonic orders.
% Taken from:
%   [1] Zotter, Frank; All-Round Ambisonic Panning and Decoding; 2012
%
% Arguments:
%   LS_pos : itaCoordinates, loudspeaker setup
%   tDesign: itaCoordinates, Ambisonics Reproduction setup
% Options: (default)
%   {N          : SH truncation order , should not be greater than floor(sqrt(L)-1)
%    'decoder'  : ('SAD')                    , Sampling Ambisonic Decoder      | Since we use a t-design this
%                 'MMAD'                     , Mode Matching Ambisonic Decoder | is actually a don't care case ([1], Appendix eq.(44))
%    'checkHull': (true)                     , check hull and calc imaginary loudspeakers
%                 false                      , do not check hull for holes
%
%   For more Options see calculate_imaginary_loudspeaker.m
%                        ita_sph_base.m , ita_vbap_pan.m
%                        spherical_harmonics.m, vbap.m
%
% Returns:
%   D_allrad: AllRAD matrix

%% Parse input
if mod(numel(varargin), 2)
    N = varargin{1};
    varargin(1) = [];
else
    N = floor(sqrt(tDesign.nPoints)-1);
end
% Default values
params.decoder = 'SAD';
params.checkHull = true;
[params, ~] = ita_parse_arguments(params, varargin);

%% Calculate VBAP matrix
% check hull for holes and fill them
if params.checkHull
    contours = ita_setup_checkForHoles(LS_pos);
    if ~isempty(contours)
        LS_imag = ita_setup_calculateImaginaryLS(LS_pos, contours);
        for idx = 1:LS_imag.nPoints
            LS_pos.cart(end+1, :) = LS_imag.cart(idx, :);
        end
    end
end

% calculate VBAP matrix
D_vbap = ita_pan_vbap(LS_pos, tDesign, 'distanceLoss', false);

%% Calculate Ambisonics Decoder
switch lower(params.decoder)
    case {'sad', 'samplingambisonicdecoder'}
        D_decoder = ita_hoa_samplingDecoder(tDesign, N, varargin{:});
    case {'mmad', 'modematchingambisonicdecoder'}
        D_decoder = ita_hoa_modeMatchingDecoder(tDesign, N, varargin{:});
    otherwise
        error('I do not know that decoder');
end

%% Calculate full Decoder Matrix
D_allrad = D_vbap * D_decoder;

% drop weights for imaginary sources
if exist('LS_imag', 'var')
    D_allrad(end - (0:LS_imag.nPoints-1), :) = [];
end

end
