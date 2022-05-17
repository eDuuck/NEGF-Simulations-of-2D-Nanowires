classdef NEGF_result
    %NEGF_RESULT A simple class to handle data transfer from calculations
    %to be analyzed in later code.
    %   NEGF_RESULT contains multiple properties that 
    %   sample points to the sample the results is related to.
    
    properties
        G {mustBeNumeric}
        Gn {mustBeNumeric}
        sigma
        sigmaIn
        sigma0 {mustBeNumeric}
        sigma0In {mustBeNumeric}
        
        reduced {mustBeNumericOrLogical}
        
        E
        sample 
    end
    
    methods
        function obj = NEGF_result(sample, E)
            %NEGF_RESULT Construct an instance of this class
            %   Detailed explanation goes here
            obj.sample = sample;
            obj.E = E;
        end
        
    end
end

