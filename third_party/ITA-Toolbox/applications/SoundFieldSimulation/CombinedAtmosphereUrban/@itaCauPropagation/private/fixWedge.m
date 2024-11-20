function [outWedge] = fixWedge(inWedge,aperturePoint)
%FIXWEDGE Makes sure wedge is located correctly


if strcmp(inWedge.edge_type,'inner_edge')
    main_face_normal = aperturePoint.main_wedge_face_normal( 1:3 );
    opposite_face_normal = aperturePoint.opposite_wedge_face_normal( 1:3 );
    aperture_start = aperturePoint.vertex_end( 1:3 );
    vertex_length = norm( aperturePoint.vertex_start - aperturePoint.vertex_end );

    outWedge = itaFiniteWedge( main_face_normal, opposite_face_normal, aperture_start, vertex_length, inWedge.edge_type );
else
    outWedge = inWedge;
end

end

