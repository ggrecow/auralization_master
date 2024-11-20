function [ data, samplerate, metadata ] = df_BRAS_directivity_ir( alpha, beta, BRAS_ir_data )
%df_BRAS_directivity_ir Note: front pole was measured on-axis of
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

if isfield( BRAS_ir_data, 'Scale' ) && isfield( BRAS_ir_data, 'NormalizeFront' )
    assert( not( BRAS_ir_data.Scale && BRAS_ir_data.NormalizeFront ) || xor( BRAS_ir_data. Scale, BRAS_ir_data.NormalizeFront ) )
end
    
% Prepare data, scale or normalize if requested
if isfield( BRAS_ir_data, 'Scale' ) && BRAS_ir_data.Scale
    
    ir = BRAS_ir_data.IR( :, idx ) ./ BRAS_ir_data.Peak;
    
elseif isfield( BRAS_ir_data, 'NormalizeFront' ) && BRAS_ir_data.NormalizeFront
    
    front = itaAudio();
    front.timeData = BRAS_ir_data.IR( :, idx );
    %front_inv = ita_invert_spk_regularization( front );
    
    ir_target = itaAudio();
    ir_target.timeData = BRAS_ir_data.IR( :, idx );
    ir_normalized = ir_target / front;
    ir = ir_normalized.timeData;
    
else
    
    ir = BRAS_ir_data.IR( :, idx ); % untouched
    
end

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

data = [ ir(1:filter_length)', zeros( append_zeros ) ];

metadata = [];
