function tf_absorbed = tf_atmospheric_absorption( obj, distance )
%tf_atmospheric_absorption Calculates atmospheric absorption after iso96713 for the distance for all frequency bins

if distance <= 0
    error 'Distance cannot be zero or negative'
end

f = obj.freq_vec( 2:end );

tf_absorbed = 1 - [ 0; ita_atmospheric_absorption_factor( f, distance ) ];  % Note: DC value set to ZERO

end
