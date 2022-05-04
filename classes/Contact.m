classdef Contact
    %CONTACT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        SC   {mustBeNumeric}
        tau {mustBeNumeric}
        beta {mustBeNumeric}
        alpha {mustBeNumeric}
        eta {mustBeNumeric}
        length
    end
    
    methods
        function obj = Contact(SC,tau)
            %CONTACT Construct an instance of this class
            %   
            if nargin == 1
                obj.SC = SC;
            elseif nargin == 2
                obj.SC = SC;
                obj.tau = tau;
            end
        end
    end
end

