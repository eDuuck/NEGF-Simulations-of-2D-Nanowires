function [remapped_data] = NEGF_result_remap(NEGF_result, data)
    %NEGF_RESULT_REMAP remaps the desired data from a NEGF_result to a matrix
    %that's the same size as the sample which can then be visualized by using
    %mehtods such as imagesc(C).
    %
    %   NEGF_RESULT_REMAP(NEGF_result) will extract and remap electron
    %   concentrations into a MxN matrix where M is the width of the sample 
    %   associated to the result and N is the length. The output can then
    %   be directly illustrated by the methdo image(C) or imagesc(C).
    %
    %   NEGF_RESULT_REMAP(NEGF_result,data) will let the desired quantity
    %   "data" be extracted from the NEGF_result. This is not case
    %   sensitive.
    %   The remapping currently supports the following quantities:
    %       "electrons"     : Just remaps Gn.
    %       "fermi"         : Extracts the fermilevels in the sample.
    %       "A"             : Remaps the spectral function which shows the
    %                       : density of states.
    %
    %
    if nargin < 2
        data = "electrons";
    end

    remapped_data = zeros(size(NEGF_result.sample.units));
    width = NEGF_result.sample.width;
    length = NEGF_result.sample.length;
    switch(lower(data))
        case "electrons"
            data_values = real(diag(NEGF_result.getGn()));
        case "fermi"
            Gn = NEGF_result.getGn();
            G = NEGF_result.getG();
            A = 1i*(G - G');
            data_values = real(diag(Gn ./A));
        case "a"
            data_values = diag(real(1i*(NEGF_result.getG() - NEGF_result.getG()')));
        otherwise
            error("Not a supported sort of data to remap.");
    end



    for j = 1:length
        remapped_data(:,j) = data_values(((j-1)*width+1):(j*width));
    end

