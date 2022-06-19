classdef NEGF_result < matlab.mixin.Copyable
    %NEGF_RESULT A simple class to handle data transfer from calculations
    %to be analyzed in later code.
    %   NEGF_RESULT contains multiple properties that 
    %   sample points to the sample the results is related to.
    
    properties
        G
        Gn
        sigma
        sigmaIn
        sigma0
        sigma0In
        
        reduced {mustBeNumericOrLogical}

        E
        B
        sample
    end

    methods
        function obj = NEGF_result(sample, E,B)
            %NEGF_RESULT Construct an instance of this class
            %   Detailed explanation goes here
            obj.sample = sample;
            obj.E = E;
            obj.B = B;
        end

        function reduce(obj,doReduce)
            if ~exist("doReduce", "var")
                doReduce = true;
            end
            if xor(doReduce, obj.reduced) %Only do calc if change is actually desired.
                if doReduce
                    obj.G = compress(obj.G,'qoi');
                    obj.Gn = compress(obj.Gn,'qoi');
                    obj.sigma0 = compress(obj.sigma0,'qoi');
                    obj.sigma0In = compress(obj.sigma0In,'qoi');
                    obj.reduced = true;
                else
                    obj.G = decompress(obj.G,'qoi');
                    obj.Gn = decompress(obj.Gn,'qoi');
                    obj.sigma0 = decompress(obj.sigma0,'qoi');
                    obj.sigma0In = decompress(obj.sigma0In,'qoi');
                    obj.reduced = false;
                end
            end
        end
    end
end

