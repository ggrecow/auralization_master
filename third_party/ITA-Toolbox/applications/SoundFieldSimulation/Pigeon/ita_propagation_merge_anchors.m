function merged_anchor = ita_propagation_merge_anchors( anchor_a, anchor_b )
%ita_propagation_merge_anchors Merges two anchors in a best-fit attempt

% Merge two outer wedges
if strcmpi( anchor_a.anchor_type, 'outer_edge_diffraction' ) && strcmpi( anchor_b.anchor_type, 'outer_edge_diffraction' )
    merged_anchor = ita_propagation_merge_outer_edge_diffraction( anchor_a, anchor_b );
else
    error( 'Merging types %s and %s is not supported, yet.', anchor_a.type, anchor_b.type )
end

end
