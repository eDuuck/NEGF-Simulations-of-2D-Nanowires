classdef Sample < handle
    %SAMPLE is a class handle that is used to describe a system that NEGF
    %simulations are to be applied to.
    %   Creating a sample class can be done by invoking 
    %   SAMPLE(w,l,eps,t,arch). w will describe the width of the wire, l
    %   describes the length, eps is the value of the basis function at a
    %   given point and t describes the interaction between different
    %   points. If t is a vector with multiple values, further neighbors
    %   are also effected. This isn't completed however and needs to be
    %   touched up upon to work correctly.
    
    properties
        width {mustBeNumeric}
        length {mustBeNumeric}
        epsilon {mustBeNumeric}
        t {mustBeNumeric}
        dim {mustBeNumeric}
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
        function obj = Sample(w,l,eps,t,arch)
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
            %up in an uncompressed raw state.
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
        
        function append(obj,M, side)
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
            obj.nbrOfContacts = obj.nbrOfContacts + 1;
            if isa(M,'Contact')
                obj.contacts{obj.nbrOfContacts} = M;
            else
                %Creates a simple 2D contact to pos.
                if nargin < 5
                    obj.contacts{obj.nbrOfContacts} = Contact(M,tau,pos);
                else
                    obj.contacts{obj.nbrOfContacts} = Contact(M,tau,pos,Inf,con_mat);
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

