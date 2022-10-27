classdef NEGF_param < matlab.mixin.Copyable
    %NEGF_PARAM Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        sample
        E
        B
        
        error_halt
        print
        errorMarg
        it_lim
        con_it_lim
        rate
        compress
        save_time 
        result

        g0
    end
    
    methods
        function obj = NEGF_param(sample,E,B,compress)
            %NEGF_PARAM Construct an instance of this class
            %   Detailed explanation goes here
            obj.sample = copy(sample);
            obj.E = E;
            if nargin < 3
                obj.B = 0;
            else
                obj.B = B;
            end
            if nargin < 4
                obj.compress = false;
            else
                obj.compress = compress;
            end
            obj.errorMarg = 1e-5;   %Seems to be accurate enough in most testing.
            obj.rate = 0.9;
            obj.it_lim = 100;
            obj.save_time = Inf;
            obj.print = false;
            obj.result = 0;
            obj.g0 = 0;     %This guess is an old NEGF_result that we expect to be similar. (Small E/B change)
            obj.error_halt = true;
            obj.con_it_lim = 10000;
        end

    end
end

