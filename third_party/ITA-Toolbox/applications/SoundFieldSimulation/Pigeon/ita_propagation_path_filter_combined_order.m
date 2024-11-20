function pps_filtered = ita_propagation_path_filter_combined_order( pps, combined_order, filter_modus )
%ita_propagation_path_length Returns paths that are within the given
%combined order (reflection + diffraction)
%   By default this does a less-equal ('leq') check. It is also possible to
%   filter using a greater-equal ('geq') check or finding a matching order
%   ('match').

if ~isfield( pps, 'propagation_anchors' ) % not a list but only one path
    error( 'Need a propagation path or path list' )
end

if nargin < 3
    filter_modus = 'leq';
end

N = numel( pps );
pps_filtered = [];

for n = 1:N

    propagation_path = pps( n );
    
    [ path_refl_order, path_diffr_order ] = ita_propagation_path_orders( propagation_path );
    path_combined_order = path_refl_order + path_diffr_order;
    
    if strcmpi( filter_modus , 'match' ) % exact match
        
        if combined_order == path_combined_order
            pps_filtered = [ pps_filtered, propagation_path ];
        end
        
    elseif  strcmpi( filter_modus , 'leq' ) % lower-equal
        
        if path_combined_order <= combined_order
            pps_filtered = [ pps_filtered, propagation_path ];
        end
        
    elseif  strcmpi( filter_modus , 'geq' ) % greater-equal
        
        if path_combined_order >= combined_order
            pps_filtered = [ pps_filtered, propagation_path ];
        end
        
    else
        error( 'Could not understand filter modus "%s", use "match", "leq" or "geq"', filter_modus )
    end
end

end
