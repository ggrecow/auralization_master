classdef AtmosphericRay < handle
%AtmosphericRay Representing an atmospheric sound propagation path (curved ray)
%
% This class represents a sound ray moving through the atmosphere. Its
% main properties are r - an itaCoordinates representing the ray path -
% and n - an itaCoordinates representing the wave normal at each point
% of r.
%
% Before start tracking this ray it is necessary to initialize the ray.
% Call init(<source vector>, <first normal vector>).
%
% Since this is a handle, caution is advised when "copying" a ray.
% Changing the "copy" will also change the original because the ray is
% represented by a handle (pointer).

% -------------------------------------------------------------------------
%                        ____  __________  _______
%                       //  / //__   ___/ //  _   |
%                      //  /    //  /    //  /_|  |
%                     //  /    //  /    //  ___   |
%                    //__/    //__/    //__/   |__|
%
% -------------------------------------------------------------------------
%                  ARTMatlab - ITAGeometricalAcoustics
%        (c) Copyright Institute of Technical Acoustics (ITA)
%  This file is part of the ARTMatlab application. Some rights reserved.
% You can find the license in the LICENSE.md file in the ARTMatlab folder.
%--------------------------------------------------------------------------
    
    %% Internal Properties
    properties(Access = private, Hidden = true)
        mPath = itaCoordinates;     %See r
        mNormals = itaCoordinates;  %See n
        mTimeStamps = [];           %See t
        
        mReflectionR;      	%Position where reflections occured
        mReflectionN;    	%Wavefront normals at the times of reflection
        mReflectionT;     	%Times of reflection
        mIdxReflection;     %Indices just before a reflection occured
        
        mSpreadingLoss;     %Spreading loss factor at end point (receiver)
    end
    
    %% Access Properties
    properties(GetAccess = public, SetAccess = private)
        isEigenray;           %Indicating whether ray is result of an eigenray search
        receiverSphereHit;    %Indicating whether eigenray search reached desired accuracy
        rayZoomingIterations; %Number of iterations used in ray zooming process during eigenray search)
    end
    
    properties(GetAccess = public, SetAccess = private, Dependent = true)
        r0;             %Initial position of the ray (itaCoordinates)
        n0;             %Initial wavefront normal (itaCoordinates)
        phi;            %Initial azimuth
        theta;          %Initial elevation
        r;              %Contains all points along the ray (itaCoordinates)
        n;              %Contains all wavefront normals along the ray (itaCoordinates)
        t;              %Vector with timestamps for each entry in r
        numPoints;      %Stored number of points
        numReflections; %Number of reflections that occured for this ray
        
        spreadingLoss;  %Spreading loss factor at ray end point (receiver)
    end
    
    properties(GetAccess = public, SetAccess = private, Hidden = true, Dependent = true)
        idxReflection;
    end
    
    %% Constructor & Initialization
    methods
        function this = AtmosphericRay(varargin)
            %If there are input arguments init() is called.
            %See init() for information on possible inputs.
            
            if(nargin == 0)
                return;
            end
            this.init(varargin{:});
        end
        
        function init(this, varargin)
            %Initialize this ray.
            %
            %Possible inputs:
            %r0:    Source position / initial position of ray (1x3-vec)
            %n0:    Initial wavefront normal (1x3-vec)
            %or
            %r0:    Source position / initial position of ray (1x3-vec)
            %theta: Initial elevation of wavefront normal (scalar)
            %phi:   Initial azimuth of wavefront normal (scalar)
            %
            if(nargin < 3 || nargin > 4)
                error([mfilename('class') ' number of input arguments must be between 2 and 3']);
            end
            
            r0 = varargin{1};
            [ok, r0] = this.limitInputToCoordinates(r0);
            if(~ok)
                error([mfilename('class') ' first input argument must be a single valued itaCoordinates or a 1x3-vector']);
            end
            
            if(nargin == 3)
                n0 = varargin{2};
                [ok, n0] = this.limitInputToCoordinates(n0);
                if(~ok)
                    error([mfilename('class') ' second input argument must be a single valued itaCoordinates or a 1x3-vector when using two input arguments']);
                end
                
            elseif(nargin == 4)
                theta = varargin{2};
                phi = varargin{3};
                
                if(~isscalar(phi) || ~isscalar(theta) || ~isnumeric(phi) || ~isnumeric(theta))
                    error([mfilename('class') ' when using three arguments the second and third argument must be angles (scalar and numeric values)']);
                end
                
                n0 = itaCoordinates([1 deg2rad(theta) deg2rad(phi)], 'sph');
                n0 = n0.sph2cart;
            end
            
            this.mPath = r0;
            this.mNormals = n0;
            this.mTimeStamps = 0;
            
            this.mReflectionR = itaCoordinates;
            this.mReflectionN = itaCoordinates;
        end
    end
    
    methods(Access = private, Static = true)
        function [ok, vec] = limitInputToCoordinates(vec)            
            ok=true;
            if(~isa(vec, 'itaCoordinates'))
                if( ~isnumeric(vec) ||  ~isequal(size(vec),[1 3]) )
                    ok=false;
                else
                    vec = itaCoordinates(vec);
                end
            elseif(vec.nPoints ~= 1)
                ok = false;
            end
        end
    end
    
      %% Load & Save
    methods
        function store(this, filename)
            % Stores this StratifiedAtmosphere object into a file using the JSON format
            fID = fopen(filename, 'w');
            try
                fwrite(fID, this.toJSON);
            catch err
                fclose(fID);
                rethrow(err);
            end
            fclose(fID);
        end
        %Creates a JSON string from this class
        jsonStr = toJSON(this)
    end
    methods(Static = true)
        function ray = load(filename)
            % Loads a StratifiedAtmosphere from file using the JSON format
            jsonStr = fileread(filename);
            ray = AtmosphericRay.parseJSON(jsonStr);
        end
        %Parses a json string to one or multpile AtmosphericRays
        ray = parseJSON(jsonStr);
    end
    
    %% Get functions
    methods
        function r0 = get.r0(this)
            if(isempty(this.mPath))
                r0 = [];
            else
                r0 = this.mPath.n(1);
            end
        end
        function n0 = get.n0(this)
            if(isempty(this.mNormals))
                n0 = [];
            else
                n0 = this.mNormals.n(1);
            end
        end
        function phi = get.phi(this)
            if(isempty(this.mNormals))
                phi = [];
            else
                phi = this.n0.phi_deg;
            end
        end
        function theta = get.theta(this)
            if(isempty(this.mNormals))
                theta = [];
            else
                theta = this.n0.theta_deg;
            end
        end
        function r = get.r(this)
            r = this.mPath;
        end
        function n = get.n(this)
            n = this.mNormals;
        end
        function t = get.t(this)
            t = this.mTimeStamps;
        end
        
        function nPoints = get.numPoints(this)
            nPoints = this.mPath.nPoints;
        end
        function nRefl = get.numReflections(this)
            nRefl = numel(this.mIdxReflection);
        end
        
        function out = get.idxReflection(this)
            out = this.mIdxReflection;
        end
        
        function out = get.spreadingLoss(this)
            out = this.mSpreadingLoss;
        end
    end
    
    %% Dependent Get Functions
    methods
        function [idxStart, idxEnd] = indicesOfReflectionOrder(this, reflOrder)
            %Returns the indices where the given reflection order starts
            %and ends. If the reflection order does not occur, both indices
            %are empty.

            if(this.numReflections < reflOrder)
                idxStart = [];
                idxEnd = [];
                return;
            end
            
            if(reflOrder == 0)
                idxStart = 1;
            else
                idxStart = this.mIdxReflection(reflOrder);
            end
            
            if(this.numReflections == reflOrder)
                idxEnd = this.numPoints;
            else
                idxEnd = this.mIdxReflection(reflOrder+1)-1;
            end
        end
        
        function [r, n, t] = reflectionParameters(this, reflOrder)
            if(reflOrder > this.numReflections)
                r = [];
                n = [];
                t = [];
                return;
            end
            
            r = this.mReflectionR.n(reflOrder);
            n = this.mReflectionN.n(reflOrder);
            t = this.mReflectionT(reflOrder);
        end
    end
    
    %% Adding Data
    methods
        function addData(this, newPoint, newNormal, newT)
            %Adds a new set of data to the ray.
            %The dataset contains of a point (1x3-vec), a wavefront normal
            %(1x3-vec) and a timestamp.
            
            idxNew = this.numPoints+1;
            
            this.mPath.cart(idxNew, :) = newPoint;
            this.mNormals.cart(idxNew, :) = newNormal;
            this.mTimeStamps(idxNew) = newT;
        end
        
        function addReflection(this, rGround, nGround, tGround)
            %Adds the relevant data for a reflection to this ray.
            %The necessary data contains a position (1x3-vec), a wavefront
            %normal (1x3-vec) and a timestamp.
            
            this.addData(rGround, nGround,tGround);
            
            nRefl = this.numReflections+1;
            
            this.mReflectionR.cart(nRefl, :) = rGround;
            this.mReflectionN.cart(nRefl, :) = nGround;
            this.mReflectionT(nRefl) = tGround;            
            this.mIdxReflection(nRefl) = this.numPoints;
        end
    end
    
    %% Public
    methods
        function path = pathLength(this, idxEnd)
            %Returns an approximation of the path length of this ray by
            %integrating over the dr/dt.
            %
            %For this purpose the ray is approximated as multiple segments,
            %defined by consecutive points.
            if nargin == 1; idxEnd = this.numPoints; end
            assert(all(idxEnd >= 1) && all(idxEnd <= this.numPoints), 'idxEnd is out of bounds.')
            
            path = zeros(size(idxEnd));
            dr = vecnorm(diff(this.r.cart), 2, 2);
            for idx=1:numel(idxEnd)
                path(idx) = sum( dr(1:(idxEnd(idx)-1)) );
            end
        end
        
        function b = raysHaveSameInitDirection(this,otherRay)
            %Returns true if this and a given ray have the same initial
            %angles theta and phi.
            b = ( abs(this.theta - otherRay.theta) < 100*eps )...
                && ( abs(this.phi - otherRay.phi) < 100*eps );
            %Hint: the comparision with 100*eps is done because numerical
            %issues occured when the angular distance between rays was very
            %small. In these cases, b became false although the rays had
            %the same direction
        end
        
        function varargout = interpolateToTime(this, timeVector, method)
            %Interpolates the data of this ray to the given timestamps.
            %Returns the sampling points r and optionally the wavefront
            %normals n.
            %   Also allows to specify the interpolation method (see
            %   interp1 for more information)
            if nargin == 2
                method = 'nearest';
            end
            
            if nargout
                varargout{1} = itaCoordinates( interp1(this.t, this.r.cart, timeVector, method, 'extrap') );
            end
            if nargout == 2
                varargout{2} = itaCoordinates( interp1(this.t, this.n.cart, timeVector, method, 'extrap') );
            end
        end
        
        function varargout = interpolateToRay(this, otherRay, method)
            %Interpolates the data of this ray to the timestamps of the
            %given ray. Returns the sampling points r and optionally the
            %wavefront normals n.
            %   Also allows to specify the interpolation method (see
            %   interp1 for more information)
            if nargin == 2
                method = 'nearest';
            end
            
            if nargout
                varargout{1} = itaCoordinates( interp1(this.t, this.r.cart, otherRay.t, method, 'extrap') );
            end
            if nargout == 2
                varargout{2} = itaCoordinates( interp1(this.t, this.n.cart, otherRay.t, method, 'extrap') );
            end
        end
        
        function reflOrder = reflectionOrderOfIdxT(this, idxT)
            % Returns the reflection order of a ray position of a given
            % index. Returns [] if the index is out of bounds.
            
            if(idxT < 1 || idxT > this.numPoints)
                reflOrder = [];
                return;
            end
            idxRefl = this.mIdxReflection;
            
            for rOrder = 1:length(idxRefl)
               if idxT<= idxRefl(rOrder)
                   reflOrder = rOrder-1;
                   return;
               end                                 
            end
            reflOrder = length(idxRefl);
        end
        
        function [idxT, err]  = closestIdxToTime(this, tCheck)
            %Returns the closest index to a given time stamp
            %
            %Input:
            %tCheck:    Time stamp to check (sec)
            %
            %Outpus:
            %idxT:      Index which is closest to the time
            %err:       The absolute error of t at idxT compared to the
            %           checked time stamp (err >= 0)
            
            [err, idxT] = min(abs(this.t-tCheck));
        end
        
        function r = pointAtTime(this, tCheck)
            %Returns the point of the ray at a given time stamp.
            %
            %If the time stamp is not stored within the time vector t, the
            %point is calculated using linear interpolation.
            %
            %Input:
            %tCheck:    Time stamp to check
            %
            %Output:
            %r:         Position at given time stamp (1x3 vector)
            
            [idxT1, absTError] = this.closestIdxToTime(tCheck);
            if(absTError == 0)
                r = this.r.cart(idxT1,:);
                return
            end
            
            if(this.t(idxT1) - tCheck > 0)
                idxT2 = idxT1-1;
            else
                idxT2 = idxT1+1;
            end
            dt = abs(this.t(idxT2) - this.t(idxT1));
            pathPortion = absTError/dt;
            
            r1 = this.r.cart(idxT1, :);
            r2 = this.r.cart(idxT2, :);            
            dr = r2-r1;
            
            r = r1 + pathPortion*dr;
        end
    end
end