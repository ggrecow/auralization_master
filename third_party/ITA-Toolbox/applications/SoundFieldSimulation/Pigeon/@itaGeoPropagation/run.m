function [ freq_data_linear, valid_paths ] = run( obj, with_progress_bar )
%RUN Calculates the transfer function (tf) of the superimposed (geometrical) propagation path in frequency domain

freq_data_linear = zeros( obj.num_bins, 1 );

if nargin < 2 || ~with_progress_bar
    with_progress_bar = false;
else
    h = waitbar( 0, 'Hold on, running propagation modeling' );
end

% Iterate over propagation paths, calculate transfer function and sum up
N = numel( obj.pps );
valid_paths = zeros( N, 1 );
for n = 1:N

    if with_progress_bar
        waitbar( n / N )
    end
    
    pp = obj.pps( n );
    [ pp_tf, valid ] = obj.tf( pp );
    if valid
        freq_data_linear = freq_data_linear + pp_tf;
        if any( isnan( pp_tf ) )
            x = obj.tf( pp );
        end
    end
    valid_paths( n ) = valid;
            
end

if with_progress_bar 
    close( h )
end

end
