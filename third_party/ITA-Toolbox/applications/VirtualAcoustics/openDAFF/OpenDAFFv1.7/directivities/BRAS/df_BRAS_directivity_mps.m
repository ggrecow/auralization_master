function [ freqs, data_complex, metadata ] = df_BRAS_directivity_mps( alpha, beta, BRAS_mps_data )
%df_BRAS_directivity_ir Note: front pole was measured on-axis of
%loudspeaker, so rotation of final DAFF file must be applied (pitch minus 90
%degree)

assert( isfield( BRAS_mps_data, 'MPS' ) )
assert( isfield( BRAS_mps_data, 'Frequency' ) )
assert( isfield( BRAS_mps_data, 'Phi' ) )
assert( isfield( BRAS_mps_data, 'Theta' ) )

if isfield( BRAS_mps_data, 'samplerate' )
    samplerate = BRAS_mps_data.samplerate;    
else
    samplerate = 44100;
end

phi = round( alpha );
theta = round( beta );

if theta == 0
    idx = 1;
elseif theta == 180
    idx = size( BRAS_mps_data.MPS, 2 );
else
    idx = and( BRAS_mps_data.Phi == phi, BRAS_mps_data.Theta == theta );
end

data_complex = BRAS_mps_data.MPS( :, idx )';
freqs = BRAS_mps_data.Frequency;

metadata = [];
