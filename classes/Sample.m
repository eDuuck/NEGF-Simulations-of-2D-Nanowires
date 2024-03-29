classdef Sample < matlab.mixin.Copyable
    %SAMPLE is a class handle that is used to describe a system that NEGF
    %simulations are to be applied to.
    %   Creating a sample class can be done by invoking 
    %   SAMPLE(w,l,eps,t,arch). w will describe the width of the wire, l
    %   describes the length, eps is the value of the basis functions at a
    %   given point and t describes the interaction between different
    %   points. If t is a vector with multiple values, further neighbors
    %   are also effected. This isn't completed however and needs to be
    %   touched up upon to work correctly as not all functions in the NEGF
    %   library features non-nearest neighbor interaction.
    %
    %   SAMPLE has a few functions to ease the construction of the device.
    %   These are:
    %       SAMPLE.append(M,side) which appends a matrix M to the channel
    %       of the device at side specified by 'side'. Side can have the
    %       values, 'u','d','l','r'. Standard value is 'r' if side isn't
    %       specified.
    %
    %       SAMPLE.applyNoise(amp,corLength) adds a noisefloor to the
    %       channel with a relative amplitude amp in relation to the max
    %       value in the channels eps value. The correlation length
    %       'corLength' is defined in amount of lattice points.
    %
    %       SAMPLE.addContact(M, tau, pos) adds a contact at the which
    %       upper left corner is at position pos. M is the epsilon values
    %       for the unit cell in the contact and needs to be a matrix with
    %       the witdth of the contact. Tau is the interacction between
    %       lattice points inbetween lattice points in the contact.
    
    properties
        width {mustBeNumeric}
        length {mustBeNumeric}
        epsilon {mustBeNumeric}
        t {mustBeNumeric}
        dim {mustBeNumeric}
        a {mustBeNumeric}
        M {mustBeNumeric}
        units {mustBeNumeric}
        D {mustBeNumeric}
        arch
        method
        compressed {mustBeNumericOrLogical}
        uniform {mustBeNumericOrLogical}
        contacts
        nbrOfContacts
    end
    
    methods
        function obj = Sample(w,l,eps,t,a,arch)
            %SAMPLE(w,l,eps,t,arch) Construct an instance of this class
            %   Creates a sample class which describes a sample with width
            %   w, length l, basis function value of eps and t. If t is a
            %   vector then the first value represents connection to
            %   nearest neighbor and second value represents second
            %   neighbor.
            if nargin >= 4
                obj.width = w;
                obj.length = l;
                obj.M = w*l;
                obj.epsilon = eps;
                obj.t = t;
                obj.units = ones(w,l)*eps;
                obj.uniform = true;
            end
            if nargin > 4
                obj.a = a;
            else
                obj.a = 10e-10;
            end
            if nargin > 5
                obj.arch = arch;
            else
                obj.arch = 'rectangular';
            end
            
            if nargin == 2
                obj.units = w;
                obj.dim = size(w);
                obj.length = obj.dim(2);
                obj.width = obj.dim(1);
                obj.M = numel(w);
                obj.epsilon = w;
                obj.t = l;
                if(min(obj.units == obj.units(1,1),[],'all') == 1)
                    obj.uniform = true;
                else
                    obj.uniform = false;
                end
                obj.arch = 'square';
            end
            obj.compressed = false;
            obj.dim = [obj.width, obj.length];
            obj.nbrOfContacts = 0;
            obj.contacts = {};
            obj.D = 0;
        end
        
        function status = compress(obj, doComp)
            %COMPRESS(doComp) Compresses the matrix that describes the
            %units that make up the sample or uncompresses if the sample is
            %allready in an compressed state if doComp is left empty.
            %Setting doComp to true guarantees that the sample gets
            %compressed while setting it to false guarantees that it ends
            %up in an uncompressed raw state. (Don't use this)
            if nargin > 1
                change = xor(doComp,obj.compressed);
            else
                change = true;
            end
            if change
                if ~obj.compressed
                    [obj.units, obj.method, obj.compressed] = compress(obj.units);
                else
                    obj.units = decompress(obj.units, obj.method);
                    obj.compressed = false;
                end
            end
            status = obj.compressed;
        end
        
        function append(obj,M,side)
            %APPEND(M,side) Appends the matrix M to the specified side of
            %the matrix that describes the sample. If side is left empty M
            %is appended to the right of the sample. Side can be set to:
            %
            %   'r' = right
            %   'u' = up
            %   'l' = left
            %   'd' = down
            if nargin < 3
                side = 'r';
            end
            obj.compress(false);
            if side == 'u' || side == 'd'
                if size(M,2) ~= obj.length
                    error('Length of M does not match sample length.')
                end
                if side == 'u'
                    obj.units = [M; obj.units];
                else
                    obj.units = [obj.units; M];
                end
                obj.width = obj.width + size(M,1);
            elseif side == 'l' || side == 'r'
                if size(M,1) ~= obj.width
                    error('Width of M does not match sample width.')
                end
                if side == 'l'
                    obj.units = [M, obj.units];
                else
                    obj.units = [obj.units, M];
                end
                obj.length = obj.length + size(M,2);
            else
                error('Side has to be either r,u,l,d');
            end
            obj.M = obj.M + numel(M);
            obj.dim = [obj.width, obj.length];
        end
        
        function applyNoise(obj,amp,corLength)
            %applyNoise(amp,corLength) adds a noisefloor to the
            %channel with a relative amplitude amp in relation to the max
            %value in the channels eps value. The correlation length
            %'corLength' is defined in amount of lattice points.
            if obj.compressed
                obj.units = decompress(obj.units, obj.method);
                obj.compressed = false;
            end
            relAmp = amp*max(abs(obj.getUnits),[],"all");
            if length(corLength) < 2
                noiseMat = rsgeng2D(max(obj.dim),max(obj.dim),relAmp,corLength);
            else
                noiseMat = rsgeng2D(max(obj.dim),max(obj.dim),relAmp,corLength(1),corLength(2));
            end
            obj.units = obj.units + noiseMat(1:obj.dim(1),1:obj.dim(2));
        end
        
        function M = getUnits(obj)
            %GETUNITS Returns the matrix M that represents the sample in
            %it's uncompressed form.
            if obj.compressed
                M = decompress(obj.units, obj.method);
            else
                M = obj.units;
            end
        end
        
        function addContact(obj, M, tau, pos, con_mat, fermi_level)
            %SAMPLE.addContact(M, tau, pos) adds a contact at the which
            %upper left corner is at position pos. M is the epsilon values
            %for the unit cell in the contact and needs to be a matrix with
            %the witdth of the contact. Tau is the interacction between
            %lattice points inbetween lattice points in the contact.
    
            obj.nbrOfContacts = obj.nbrOfContacts + 1;
            if isa(M,'Contact')
                obj.contacts{obj.nbrOfContacts} = M;
            else
                %Creates a simple 2D contact to pos.
                if nargin < 5
                    obj.contacts{obj.nbrOfContacts} = Contact(M,tau,pos,obj.a);
                else
                    obj.contacts{obj.nbrOfContacts} = Contact(M,tau,pos,obj.a,Inf,con_mat);
                end
                if pos(2) == 1
                    obj.contacts{obj.nbrOfContacts}.face = 1;
                else 
                    obj.contacts{obj.nbrOfContacts}.face = -1;
                end
                if nargin < 6
                    obj.contacts{obj.nbrOfContacts}.fermi = 1;
                else
                    obj.contacts{obj.nbrOfContacts}.fermi = fermi_level;
                end
            end
        end
    end
end

