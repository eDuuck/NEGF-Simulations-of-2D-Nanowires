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
    %       "electrons"     : Just remaps Gn (divided by 2pi).
    %       "fermi"         : Extracts the fermilevels in the sample.
    %       "A"             : Remaps the spectral function which shows the
    %                       : density of states.
    %
    %
    if nargin < 2
        data = "electrons";
    end

    remapped_data = zeros(size(NEGF_result.sample.units));
    wid = NEGF_result.sample.width;
    len = NEGF_result.sample.length;
    switch(lower(data))
        case "electrons"
            data_values = real(diag(NEGF_result.getGn()))/(2*pi);
        case "fermi"
            G = NEGF_result.getG();
            sigIn = NEGF_result.getSigmaIn();
            sigInSum = zeros(size(G));
            for j = 1:length(sigIn)
                sigInSum = sigInSum + sigIn{j};
            end
            Gn = G *(sigInSum + NEGF_result.getSigma0In()) * G';

            A = 1i*(G - G');
            data_values = real(diag(Gn ./A));
        case "a"
            data_values = diag(real(1i*(NEGF_result.getG() - NEGF_result.getG()')))/(2*pi);
        otherwise
            error("Not a supported sort of data to remap.");
    end

    for j = 1:len
        remapped_data(:,j) = data_values(((j-1)*wid+1):(j*wid));
    end

