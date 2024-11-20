function [filterMags,filterGain] = calcDiffractionFilter(urbanPath, freqVector, c, maxReflOrder)
%CALCDIFFRACTIONFILTER Calculates diffraction filter for given urban path

filterMags = ones( size(freqVector) );
filterGain = 1;

propagation_anchors = urbanPath.propagation_anchors;
nAnchors = numel(propagation_anchors);

if nAnchors == 2
    return
end

number_of_diff = 0;

% check number of reflection is max reflection order
reflection_order = ita_propagation_path_orders( urbanPath );
max_refl_reached = reflection_order >= maxReflOrder;

for i = 2:nAnchors-1
    
    a_prev = propagation_anchors{i-1};
    a_curr = propagation_anchors{i};
    a_next = propagation_anchors{i+1};
    
    anchor_type = a_curr.anchor_type;
    
    segment_distance = norm( a_curr.interaction_point - a_prev.interaction_point );
    
    if strcmp(anchor_type,'outer_edge_diffraction') || strcmp(anchor_type,'inner_edge_diffraction') 
              
        source_pos = a_prev.interaction_point(1:3);
        target_pos = a_next.interaction_point(1:3);

        w = ita_propagation_wedge_from_diffraction_anchor( a_curr ); 
        if ~w.point_on_aperture(a_curr.interaction_point(1:3))
            w = fixWedge(w,a_curr);
        end       
        w.set_get_geo_eps( 1e-6 );
        w.location = a_curr.interaction_point(1:3);
%         w.set_get_geo_eps( 1e-4 );
%         precision = abs(log10(w.set_get_geo_eps( 1e-4 )))-1;
%         w.location = round(a_curr.interaction_point(1:3),precision);
                
        rho = ita_propagation_effective_source_distance( urbanPath, i ); %effective distance from aperture point to source
        last_pos_dirn = a_prev.interaction_point(1:3) - a_curr.interaction_point(1:3); %direction to the last source
        eff_source_pos = ( last_pos_dirn .* rho ./ norm(last_pos_dirn) ) + a_curr.interaction_point(1:3);
        r = ita_propagation_effective_target_distance( urbanPath, i ); %effective distance from aperture point to receiver
        next_pos_dirn = a_next.interaction_point(1:3) - a_curr.interaction_point(1:3); %"receiver"
        eff_receiver_pos = ( next_pos_dirn .* r ./ norm(next_pos_dirn) ) + a_curr.interaction_point(1:3);

        aperture_point = w.get_aperture_point( source_pos, target_pos );

        try 
            [~, D, A] = ita_diffraction_utd( w, eff_source_pos, eff_receiver_pos, freqVector, c, aperture_point );
            
            assert( ~isnan(D(1)) )
            assert( D(1) ~= 0 )
            % set correct sign
            if eff_source_pos(3) < 0 %|| eff_receiver_pos(3) < 0
                use_main_face = true;
            else
                use_main_face = false;
            end
            filterGain = filterGain*getPhaseSign(w,eff_source_pos,eff_receiver_pos,use_main_face,max_refl_reached);
        catch        
            % correct sign an diffraction filter could not be determined
            % -> neglect path and assume positive sign
            D = 1e-6*ones(size(freqVector));
            A = 1e-6;
        end                      
               
        number_of_diff = number_of_diff + 1;
        if( number_of_diff == 1 )
            filterGain = filterGain * (A / rho);
        else
            filterGain = filterGain * A;
        end
        filterMags = filterMags .* D;   
        
    else
        rho = ita_propagation_effective_source_distance( urbanPath, i ); 
        filterGain = filterGain *1/segment_distance;
        
        % if current segment is last segment also consider distance to
        % receiveiver
        if strcmp(a_next.anchor_type,'receiver')
            filterGain = filterGain*1/norm( a_next.interaction_point - a_curr.interaction_point );
        end
    end
    
end
