classdef NEGF_param < matlab.mixin.Copyable
    %NEGF_PARAM is a data object used to specify the parameters in NEGF.m,
    %NEGF_map.m and sigma_from_sample.m. A specified sample has to be fed
    %into NEGF_param along with a energy value or range E. Other parameters
    %has standard values.
    %   NEGF_PARAM(sample,E) will create a NEGF_PARAM object with standarad
    %   parameters connected to the sample and with an energy range
    %   specified by E.
    %
    %   NEGF_PARAM(sample,E,B) will also specify a magnetic field value
    %   range according to B. If B is left out then B will be set to 0.
    %
    %   NEGF_PARAM(sample,E,B,compress) will let the user specify if the
    %   result should be compressed after performing the NEGF calculations.
    %   'compress' should be set to 'true' or 'false'. 
    %
    %   The deafult values for all parameters are:
    %       sample          -       Undefined
    %       E               -       Undefined
    %       B               -       0
    %    
    %       error_halt      -       True
    %       print           -       False
    %       errorMarg       -       1E-5
    %       it_lim          -       100
    %       con_it_lim      -       1000
    %       rate            -       0.9
    %       compress        -       False
    %       save_time       -       Inf
    %
    %       g0              -       0
    %
    %   These values are used for:
    %       errorMarg   : The margin of error that is allowed from the
    %                     converging calculations.
    %
    %       it_lim      : Maximum allowed iterations for the converging
    %                     calculations used to calculate Sigma_0 in NEGF.m.
    %
    %       con_it_lim  : Maximum allowed iterations for the converging
    %                     calculations used to calculate Sigma in
    %                     contact_surface.
    %
    %       error_halt  : If the algorithms halt if the converging 
    %                     calculations did not converge according to
    %                     errorMarg, it_lim and con_it_lim.
    %   
    %       rate        : The rate of change in the converging algorithms.
    %
    %       save_time   : How long the calculations are run before the
    %                     NEGF_map.m function saves the current process of
    %                     the simulation.
    %       
    %       compress    : If the results from NEGF.m should be compressed
    %                     before being returned.
    %
    %       print       : If the progress should be printed in the command
    %                     window when using NEGF_map.m.
    %
    %       g0          : A guess of what the result should be in NEGF.
    %                     This is used in NEGF_map.m to improve calculation
    %                     times.
    
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

        result %This is used by NEGF_map to continue autosaved results.

        g0
    end
    
    methods
        function obj = NEGF_param(sample,E,B,compress)
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
            obj.result = 0; %Honestly I have no idea what this is used for.
            obj.g0 = 0;     %This guess is an old NEGF_result that we expect to be similar. (Small E/B change)
            obj.error_halt = true;
            obj.con_it_lim = 1000;
        end

    end
end

