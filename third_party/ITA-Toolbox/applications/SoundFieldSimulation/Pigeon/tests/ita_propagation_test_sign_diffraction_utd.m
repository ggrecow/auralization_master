%
% When auralising with energetic filters, the sign can be used to invert
% the phases of a path contributing to the overall field. General rule:
%
% - Direct sound paths and reflected sound paths are not flipping the sign
%   (i.e., have positive sign "+1")
% - Diffracted sound changes phase after crossing the illumination
%   boundaries
%
% This requires a "sign barrier" between the two boundaries, where the sign
% is flipped. One option is to determine the barrier geometrically "in the
% middle" of the shadow and reflection boundary.
%
% Because the problem is symmetrical, source and receiver can be exchanged.
%

n1 = [ 0 -1 0 ];
n2 = [ -1 0 0 ];
loc = [ 0 0 0 ];

w1 = itaInfiniteWedge( n1, n2, loc );
w2 = itaInfiniteWedge( n2, n1, loc );

R = [ -1 0.9 0 ];
S = [ 1 -1 0 ];

assert( calc_sign( w1, R, S ) == calc_sign( w1, S, R ) )
assert( calc_sign( w1, R, S ) == calc_sign( w2, R, S ) )
assert( calc_sign( w2, R, S ) == calc_sign( w1, R, S ) )
assert( calc_sign( w2, R, S ) == calc_sign( w2, S, R ) )

function g_sign = calc_sign( w, R, S )

if ita_diffraction_shadow_zone( w, R, S )
    
    g_sign = +1.0; % if inside shadow zone, always use positive sign
    
else
    
    both_entities_on_same_side = false;
    
    barrier = w.opening_angle / 2;
    if ( w.angle_main_face( S ) <= barrier && w.angle_main_face( R ) <= barrier )
        reflection_at_main_face = true;
        both_entities_on_same_side = true;
    elseif ( w.angle_main_face( S ) > barrier && w.angle_main_face( R ) > barrier )
        reflection_at_main_face = false;
        both_entities_on_same_side = true;
    end
    
    if both_entities_on_same_side
        
        if ita_diffraction_reflection_zone( w, R, S, reflection_at_main_face )
            g_sign = -1.0;
        else
            g_sign = +1.0;
        end
        
    else % mixed situation ... complicated
        
        reflection_at_main_face = false;
        flip_entities_symmetrical_solution = false;
        if ( w.angle_main_face( S ) + w.angle_main_face( R ) ) / 2 < barrier
            reflection_at_main_face = true;
            if w.angle_main_face( S ) < w.angle_main_face( R )
                flip_entities_symmetrical_solution = true;
            end
        else
            if w.angle_main_face( S ) > w.angle_main_face( R )
                flip_entities_symmetrical_solution = true;
            end
        end

        if flip_entities_symmetrical_solution
            g_sign = ita_diffraction_utd_illumination_sign(  w, S, R, reflection_at_main_face );
        else
            g_sign = ita_diffraction_utd_illumination_sign(  w, R, S, reflection_at_main_face );
        end
        
    end
end

end
