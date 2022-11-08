classdef Contact
    %CONTACT is a simple dataclass containing the required data of a
    %contact connected to a sample.
    
    properties
        SC {mustBeNumeric}
        tau {mustBeNumeric}
        eta {mustBeNumeric}
        unit_length
        a {mustBeNumeric}
        fermi {mustBeNumeric}
        con_mat {mustBeNumeric}
        pos (1,2) {mustBeNumeric}

        face
    end
    
    methods
        function obj = Contact(SC,tau,pos,a,unit_length,con_mat)
            %CONTACT Construct an instance of this class
            %
            obj.SC = SC;
            obj.tau = tau;
            obj.fermi = 1;
            M = numel(SC);
            t_l = length(tau);
            if nargin < 5
                obj.unit_length = Inf;
            else
                obj.unit_length = unit_length;
            end
            if nargin < 6
                obj.con_mat = zeros(M*t_l);
                for j = 1:t_l
                    obj.con_mat(1:M,(M*(j-1)+1):M*j) = eye(M) * tau(1);
                end
            else
                obj.con_mat = con_mat;
            end
            obj.eta = min(abs(tau),[],'all')*1e-9;
            obj.pos = pos;
            obj.a = a;
        end
    end
end


