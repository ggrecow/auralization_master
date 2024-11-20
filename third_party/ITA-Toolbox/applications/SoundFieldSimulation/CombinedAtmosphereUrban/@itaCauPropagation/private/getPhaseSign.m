function [outSign] = getPhaseSign(wedge,source_pos,receiver_pos,use_main_face,ignore_refl_zone)
%UNTITLED determines correct sign of diffraction path depending on relative
%phase to the other paths

if nargin < 4 
    use_main_face = false;
end

if nargin < 5
    ignore_refl_zone = false;
end

outSign = 1;

if ~isequal(wedge.main_face_normal',[0 0 1]) && ~isequal(wedge.opposite_face_normal',[0 0 1])
    % tests do not work properly, assume positve sign
    return
end
%% check if inside shadow zone
if ita_diffraction_shadow_zone(wedge,source_pos,receiver_pos)
    return
end
%% check if inside reflection zone
if ita_diffraction_reflection_zone(wedge,source_pos,receiver_pos,use_main_face) 
    if ignore_refl_zone
        outSign = 0;
    else
        outSign = -1;
    end
    return
end
%% chose sign in illuminated region depedning on distances to wedges

main_face_up = isequal(wedge.main_face_normal',[0 0 1]);
outSign = ita_diffraction_utd_illumination_sign(wedge,source_pos,receiver_pos,main_face_up);
    
if outSign == 1 && ignore_refl_zone
    %ignore path
    outSign = 0;
end

end

