classdef Contact
    %CONTACT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        SC {mustBeNumeric}
        tau {mustBeNumeric}
        beta {mustBeNumeric}
        alpha {mustBeNumeric}
        eta {mustBeNumeric}
        basis_length
        fermi {mustBeNumeric}
        con_mat {mustBeNumeric}
        pos (1,2) {mustBeNumeric}
    end
    
    methods
        function obj = Contact(SC,tau,pos,basis_length,con_mat)
            %CONTACT Construct an instance of this class
            %
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
            if nargin < 4
                obj.basis_length = Inf;
            else
                obj.basis_length = basis_length;
            end
            if nargin < 5
                obj.con_mat = zeros(M*t_l);
                for j = 1:t_l
                    obj.con_mat(1:M,(M*(j-1)+1):M*j) = eye(M) * tau(1);
                end
            else
                obj.con_mat = con_mat;
            end
            obj.eta = min(abs(tau),[],'all')*1e-3;
            obj.pos = pos;
        end
    end
end


