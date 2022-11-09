classdef NEGF_result < matlab.mixin.Copyable
    %NEGF_RESULT A simple class to handle data transfer from calculations
    %to be analyzed in later code. NEGF_RESULT contains multiple properties
    %that is used for future calculations and data analysis. NEGF_RESULT 
    %also contains functions to compress the data and extract the 
    %compressed data. To retrieve the data use the NEGF_RESULT.getX
    %functions.
    %
    %   NEGF_RESULT.getG() returns the greens function matrix.
    %
    %   NEGF_RESULT.getGn() returns the electron density matrix.
    %
    %   NEGF_RESULT.getA() returns the spectral function matrix.
    %
    %   NEGF_RESULT.getSigma() returns the Sigmas in a cell array. If the
    %   Sigma matrix for a specific contact is desired, use 
    %   NEGF_RESULT.getSigma(contact).
    %
    %   NEGF_RESULT.getSigmaIn() returns the SigmaIn in a cell array. If 
    %   the SigmaIn matrix for a specific contact is desired, use 
    %   NEGF_RESULT.getSigmaIn(contact).
    %
    %   NEGF_RESULT.getSigma0() returns the scattering matrix Sigma0.
    %
    %   NEGF_RESULT.getSigma0In() returns the scattering matrix Sigma0In.
    %
    %   NEGF_RESULT.getS0(contact) returns the guess for the contact SGF.
    %
    %   NEGF_RESULT.compress() compresses some of the matrices in the
    %   NEGF result. This is a lossy compression and will currently reduce
    %   the resolution of the data to 8 bits. On top of this loss-less
    %   compression (QOI) is used to further compress data and using the
    %   symmetry of cerrtain matrices to further save memory. A reduction
    %   of at least 8x should be expected. If uncompression is desired, use
    %   NEGF_result.compress(false).

    properties
        %G
        sigma
        %sigmaIn
        sigma0
        sigma0In
        fermiLevels

        s0  %Guess for surface function.

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
            obj.sample = sample;
            obj.E = E;
            obj.B = B;
            obj.reduced = false;
            obj.compressionMethod = 'QOI';
        end

        function compress(obj,doReduce)
            %   Detailed explanation goes here
            if ~exist("doReduce", "var")
                doReduce = true;
            end
            if xor(doReduce, obj.reduced) %Only do comp if change is actually desired.
                if doReduce
                    %obj.G = compress(obj.G, obj.compressionMethod);
                    for i = 1:length(obj.sigma)
                        obj.sigma{i} = compress(obj.sigma{i}, obj.compressionMethod);
                        obj.s0{i} = compress(obj.s0{i}, obj.compressionMethod);
                    end
                    obj.sigma0 = compress(obj.sigma0, obj.compressionMethod);
                    obj.sigma0In = compress(obj.sigma0In, obj.compressionMethod);
                    obj.reduced = true;
                else
                    %obj.G = decompress(obj.G, obj.compressionMethod);
                    for i = 1:length(obj.sigma)
                        obj.sigma{i} = decompress(obj.sigma{i}, obj.compressionMethod);
                        obj.s0{i} = decompress(obj.s0{i}, obj.compressionMethod);
                    end
                    obj.sigma0 = decompress(obj.sigma0, obj.compressionMethod);
                    obj.sigma0In = decompress(obj.sigma0In, obj.compressionMethod);
                    obj.reduced = false;
                end
            end
        end


        function G = getG(obj)
            EI = obj.E*eye(obj.sample.M);
            H = hamiltonian(obj.sample,obj.B);
            sig = obj.getSigma();
            sigSum = zeros(size(sig{1}));
            for j = 1:length(sig)
                sigSum = sigSum + sig{j};
            end
            sigSum = sigSum + obj.getSigma0();

            G = (EI - H - sigSum)^-1;

%             if obj.reduced
%                 G = decompress(obj.G, obj.compressionMethod);
%             else
%                 G = obj.G;
%             end
        end

        function Gn = getGn(obj)
            Gf = obj.getG(); 
            sigIn = obj.getSigmaIn();
            sigInSum = zeros(size(Gf));
            for j = 1:length(sigIn)
                sigInSum = sigInSum + sigIn{j};
            end
            Gn = Gf *(sigInSum + obj.getSigma0In()) * Gf';
        end

        function A = getA(obj)
            Gf = obj.getG();
            A = 1i*(Gf - Gf');
        end

        function sigma = getSigma(obj,contact)
            if nargin < 2
                sigma = cell(size(obj.sigma));
                if obj.reduced
                    for i = 1:length(obj.sigma)
                        sigma{i} = decompress(obj.sigma{i},obj.compressionMethod);
                    end
                else
                    for i = 1:length(obj.sigma)
                        sigma{i} = obj.sigma{i};
                    end
                end
            else
                if obj.reduced
                    for i = 1:length(obj.sigma)
                        sigma = decompress(obj.sigma{contact},obj.compressionMethod);
                    end
                else
                    for i = 1:length(obj.sigma)
                        sigma = obj.sigma{contact};
                    end
                end
            end
        end

        function sigmaIn = getSigmaIn(obj,contact)
            if nargin < 2
               sigmaIn = obj.getSigma();
               for i = 1:length(obj.sigma)
                    sigmaIn{i} = obj.fermiLevels(i)*1i*(sigmaIn{i}-sigmaIn{i}');
               end
            else
                sig = obj.getSigma{contact};
                sigmaIn = obj.fermiLevels(contact)*1i*(sig-sig');
            end
        end

        function sigma0In = getSigma0In(obj)
            if obj.reduced
                sigma0In = decompress(obj.sigma0In, obj.compressionMethod);
            else
                sigma0In = obj.sigma0In;
            end
        end

        function sigma0 = getSigma0(obj)
            if obj.reduced
                sigma0 = decompress(obj.sigma0, obj.compressionMethod);
            else
                sigma0 = obj.sigma0;
            end
        end

        function gamma = getGamma(obj, contact)
            if ~exist("contact","var")
                sig = obj.getSigma();
                sigSum = zeros(size(sig));
                for j = 1:length(obj.femiLevels)
                    sigSum = sigSum + sig{j};
                end
                gamma = 1i*(sigSum-sigSum');
                return
            end
            sig = obj.getSigma(contact);
            gamma = 1i*(sig-sig');
        end
        
        function s0 = getS0(obj,contact)
            if nargin < 2
                s0 = cell(size(obj.s0));
                if obj.reduced
                    for i = 1:length(obj.s0)
                        s0{i} = decompress(obj.s0{i},obj.compressionMethod);
                    end
                else
                    for i = 1:length(obj.s0)
                        s0{i} = obj.s0{i};
                    end
                end
            else
                if obj.reduced
                    for i = 1:length(obj.s0)
                        s0 = decompress(obj.s0{contact},obj.compressionMethod);
                    end
                else
                    for i = 1:length(obj.s0)
                        s0 = obj.s0{contact};
                    end
                end
            end
        end
    end
end