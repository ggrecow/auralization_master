function art = initializeRayTracer()
%INITIALIZERAYTRACER 

art = AtmosphericRayTracer;  

art.abortMaxNAdaptations = 50;
art.abortMinAngleResolution = 1e-5;

art.advancedRayZoomingThreshold = 0.1;

art.maxAngleForGeomSpreading = 0.01;
art.maxPropagationTime = 50;
art.maxReceiverRadius = 0.025;    

art.maxReflectionOrder = 0;
art.maxSourceReceiverAngle = 1;

end
