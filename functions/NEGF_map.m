function output = NEGF_map(NEGF_param)
% NEGF_MAP is a more flexible function that NEGF.m, being able to handle a
% range of values for E and B and returns a structure containing
% NEGF_results for every combination.
%   The returned structure will contain the following properties:
%
%       NEGF_result :   A 2D cell with the results from the NEGF
%                       in the form of NEGF_result{E,B}.
%
%       completed   :   A 2D matrix showing if the correlating NEGF_result
%                       at the same index. A completed simulation changes
%                       the value at completed(E,B) to 1.
%
%       E           :   The values of E.
%
%       B           :   The values of B.

save_time = NEGF_param.save_time;

E = NEGF_param.E;
B = NEGF_param.B;

if NEGF_param.result ~= 0
    output = NEGF_param.result;
else
    output = struct("E",E,...
    "B",B,"completed", zeros(length(B),length(E)));
    output.NEGF_result = cell(length(B),length(E));
end

saveTimer = tic; %Used to monitor time for autosave.
if NEGF_param.print
    disp(0 + "/" + length(B)*length(E));
end
for i = 1:length(E)
    for j = 1:length(B)
        if ~output.completed(j,i)
            params = copy(NEGF_param); %This is so that when we select value for E/B, we won't change the value of NEGF_param.
            params.E = E(i); params.B = B(j);

            if j ~= 1
                params.g0 = output.NEGF_result{j-1,i};
            elseif i ~= 1 
                params.g0 = output.NEGF_result{j,i-1};
            end

            output.NEGF_result{j,i} = NEGF(params);
            output.completed(j,i) = 1;
            if toc(saveTimer) > save_time
                %saveData %This functionality hasn't been implemented yet.
                saveTimer = tic;
            end
            if NEGF_param.print
                disp(length(B)*(i-1)+j + "/" + length(B)*length(E));
            end
        end
    end
end
end

