classdef Contact
    %CONTACT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        SC   {mustBeNumeric}
        tau {mustBeNumeric}
        beta {mustBeNumeric}
        alpha {mustBeNumeric}
        eta {mustBeNumeric}
        size
    end
    
    methods
        function obj = Contact(SC,tau,size)
            %CONTACT Construct an instance of this class
            %
            if nargin == 1
                obj.SC = SC;
            elseif nargin > 1
                obj.SC = SC;
                obj.tau = tau;
                M = numel(SC);
                t_l = length(tau);
                temp = Sample(M,t_l,SC(1),tau);
                obj.alpha = hamiltonian(temp);
                obj.beta = zeros(t_l*M);
                for j = 1:t_l
                    obj.beta = obj.beta + ...
                        diag(ones(M*(t_l-j+1),1),M*(j-1))*tau(t_l-j+1);
                end
                obj.beta = sparse(obj.beta);
                obj.size = size;
                obj.eta = min(abs(tau),[],'all')*1e-3;
            end
        end
    end
end

