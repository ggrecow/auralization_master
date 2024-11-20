function bin_out=ita_hoa_binaural_test(pos,directPlay)
% Encodes and decodes a virtual sound source position to the desired position
% and decodes it using the vrlab setup. Afterwards it is mix down to a
% binaural audio object. pos has to be a 1x3 array. play indicates whether
% a direct playback is desired.

% Coresponding: MKO@akustik.rwth-aachen.de
if nargin<1
    error('need a virtual source position');
end
if nargin<2
    directPlay=true;
end

vs_signal=ita_demosound;
vs_signal.trackLength=4;
vs_signal.channelCoordinates.cart=pos;

bformat=ita_hoa_encode(vs_signal,2);
ls_signals=ita_hoa_decode(bformat,ita_setup_LS_VRLab('virtualSpeaker',true),'vrlab',true);
bin_out=ita_binauralMixdown(ls_signals,'LSPos',ita_setup_LS_VRLab);

if directPlay
    bin_out.play
end
bar(ls_signals.rms);
end

