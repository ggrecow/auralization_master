function [ freqs, data, metadata ] = df_BRAS_directivity_ms( alpha, beta, BRAS_mps_data )

[ freqs, data_complex, metadata ] = df_BRAS_directivity_mps( alpha, beta, BRAS_mps_data );
data = abs( data_complex );
