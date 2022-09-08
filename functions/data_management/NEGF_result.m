classdef NEGF_result < matlab.mixin.Copyable
    %NEGF_RESULT A simple class to handle data transfer from calculations
    %to be analyzed in later code. NEGF_RESULT contains multiple properties
    %that is used for future calculations and data analysis. NEGF_RESULT 
    %also contains functions to compress the data and extract the 
    %compressed data.
    %
    %   NEGF_RESULT.getG() returns the greens function matrix from the NEGF
    %   result and should always be used as it will decompress G if it
    %   allready is compressed.
    %
    %   NEGF_RESULT.reduce() compresses some of the matrices in the
    %   NEGF result. This is a lossy compression and will currently reduce
    %   the resolution of the data to 8 bits. On top of this loss-less
    %   compression (QOI) is used to further compress data and using the
    %   symmetry of cerrtain matrices to further save memory. A reduction
    %   of at least 8x should be expected.

    properties
        G
        sigma
        sigmaIn
        sigma0
        sigma0In
        fermiLevels

        g0

        compressionMethod

        E
        B
        sample
    end

    properties (Access = private)
        reduced {mustBeNumericOrLogical}
    end

    methods
        function obj = NEGF_result(sample, E,B)
            %NEGF_RESULT Construct an instance of this class
            %   Detailed explanation goes here
            obj.sample = sample;
            obj.E = E;
            obj.B = B;
            obj.reduced = false;
            obj.compressionMethod = 'gomp';
        end

        function reduce(obj,doReduce)
            %REDUCE Lmao, you really reading this?
            %   Detailed explanation goes here
            if ~exist("doReduce", "var")
                doReduce = true;
            end
            if xor(doReduce, obj.reduced) %Only do comp if change is actually desired.
                if doReduce
                    obj.G = compress(obj.G, obj.compressionMethod);
                    obj.sigma0 = compress(obj.sigma0, obj.compressionMethod);
                    obj.sigma0In = compress(obj.sigma0In, obj.compressionMethod);
                    obj.reduced = true;
                else
                    obj.G = decompress(obj.G, obj.compressionMethod);
                    obj.sigma0 = decompress(obj.sigma0, obj.compressionMethod);
                    obj.sigma0In = decompress(obj.sigma0In, obj.compressionMethod);
                    obj.reduced = false;
                end
            end
        end


        function G = getG(obj)
            EI = obj.E*eye(obj.sample.M);
            H = hamiltonian(obj.sample,obj.B);

            sigSum = zeros(size(obj.sigma{1}));
            for j = 1:length(obj.sigma)
                sigSum = sigSum + obj.sigma{j};
            end
            sigSum = sigSum + obj.sigma0;

            G = (EI - H - sigSum)^-1;

%             if obj.reduced
%                 G = decompress(obj.G, obj.compressionMethod);
%             else
%                 G = obj.G;
%             end
        end


        function sigmaIn = getSigmaIn(obj)
            sigmaIn = obj.sigma;
            for i = 1:obj.fermiLevels
                sigmaIn{i} = sigmaIn{i}*obj.fermiLevels(i);
            end
        end

        function Gn = getGn(obj)
            Gf = obj.getG();
            sigInSum = zeros(size(Gf));
            for j = 1:length(obj.sigmaIn)
                sigInSum = sigInSum + 1i*(obj.sigmaIn{j} - obj.sigmaIn{j}');
            end
            Gn = Gf *(sigInSum + obj.getSigma0In()) * Gf';
        end

        function A = getA(obj)
            Gf = obj.getG();
            A = 1i*(Gf - Gf');
        end

        function sigma0In = getSigma0In(obj)
            if obj.reduced
                sigma0In = decompress(obj.sigma0In, obj.compressionMethod);
            else
                sigma0In = obj.sigma0In;
            end
        end

        function gamma = getGamma(obj, contact)
            if ~exist("contact","var")
                sigSum = zeros(size(obj.sigma));
                for j = 1:length(obj.femiLevels)
                    sigSum = sigSum + obj.sigma{j};
                end
                gamma = 1i*(sigSum-sigSum');
                return
            end
            gamma = 1i*(obj.sigma{contact}-obj.sigma{contact}');
        end

    end
end