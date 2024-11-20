function [ data, samplerate, is_symetric, metadata ] = df_BRAS_directivity_dft( alpha, beta, BRAS_ir_data )
%df_BRAS_directivity_dft Note: front pole was measured on-axis of
%loudspeaker, so rotation of final DAFF file must be applied (pitch minus 90
%degree)

assert( isfield( BRAS_ir_data, 'IR' ) )
assert( isfield( BRAS_ir_data, 'Phi' ) )
assert( isfield( BRAS_ir_data, 'Theta' ) )

if isfield( BRAS_ir_data, 'samplerate' )
    samplerate = BRAS_ir_data.samplerate;    
else
    samplerate = 44100;
end

phi = round( alpha );
theta = round( beta );

if theta == 0
    idx = 1;
elseif theta == 180
    idx = size( BRAS_ir_data.IR, 2 );
else
    idx = and( BRAS_ir_data.Phi == phi, BRAS_ir_data.Theta == theta );
end

ir = BRAS_ir_data.IR( :, idx );
ir_length = numel( ir );
append_zeros = 0;

if isfield( BRAS_ir_data, 'filter_length' )
    if BRAS_ir_data.filter_length > ir_length
        append_zeros = BRAS_ir_data.filter_length - ir_length;
        filter_length = ir_length;
    else
        filter_length = BRAS_ir_data.filter_length;
    end
else
    filter_length = ir_length;
end

% DAFF requires data alignment by multiple of 4
append_zeros = append_zeros + mod( filter_length, 4 );

ir_data_padded = [ ir(1:filter_length)', zeros( append_zeros ) ];
dft_data = fft( ir_data_padded );
N = size( dft_data, 2 );
data = dft_data( 1 : N / 2 + 1 ) / N * sqrt( 2 );

is_symetric = true;
metadata = [];
