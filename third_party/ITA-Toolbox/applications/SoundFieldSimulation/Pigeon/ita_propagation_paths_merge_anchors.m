function pps_merged = ita_propagation_paths_merge_anchors( pps, merging_gap )
%ita_propagation_path_merge_anchors Returns an altered set of propagation
% pats with potentially merged anchors where distances between paths are 
% smaller than the merging gap

if ~isfield( pps, 'propagation_anchors' ) % not a list but only one path
    error( 'Need a propagation path or path list' )
end

if nargin < 2
    merging_gap = 17e-3; % wave length of 20 kHz at 340 m/s
end

assert( merging_gap > 0 )

N = numel( pps );

num_merged = 0;

for n=1:N
    
    propagation_path = pps( n );
    M = numel( propagation_path.propagation_anchors );
    
    if M == 2
        dist_vec = propagation_path.propagation_anchors( 2 ) - propagation_path.propagation_anchors( 1 );
        if norm( dist_vec ) < merging_gap
            warning( 'Direct sound path smaller than mering gap, skipped in merged paths' );
        else
            num_merged = num_merged + 1;
            pps_merged( num_merged ) = propagation_path;
        end
        continue
    end
    
    propagation_path_merged = propagation_path;
    propagation_path_merged.propagation_anchors = propagation_path.propagation_anchors( 1 );
    
    k = 0; % number of merged anchors
    m = 1;
    for m = 1:M-1
        cur_anchor = propagation_path_merged.propagation_anchors{ m - k };
        next_anchor = propagation_path.propagation_anchors{ m + 1 };
        dist_vec = next_anchor.interaction_point( 1:3 ) - cur_anchor.interaction_point( 1:3 );
        if norm( dist_vec ) < merging_gap
            % Potentially override old merged anchor
            propagation_path_merged.propagation_anchors{ m - k } = ita_propagation_merge_anchors( cur_anchor, next_anchor );
            k = k + 1;
        else
            propagation_path_merged.propagation_anchors{ m + 1 - k } = next_anchor;
        end
    end
    
    num_merged = num_merged + 1;
    pps_merged( num_merged ) = propagation_path_merged;
    
end

end
