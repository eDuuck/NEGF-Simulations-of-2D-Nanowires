function output = NEGF_map(NEGF_param)
% NEGF_map(sample,E,B,compress,result,saveTime,errorMarg,rate,it_lim)
%   Detailed explanation goes here
save_time = NEGF_param.save_time;

E = NEGF_param.E;
B = NEGF_param.B;

newRes = true;  %If an old simulation that was autosaved is available.
if NEGF_param.result ~= 0
    output = result;
    newRes = false;
end 

if newRes
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
            params = copy(NEGF_param);
            params.E = E(i); params.B = B(j);

            if j ~= 1
                params.g0 = output.NEGF_result{j-1,i};
            elseif i ~= 1 
                params.g0 = output.NEGF_result{j,i-1};
            end

            output.NEGF_result{j,i} = NEGF(params);
            output.completed(j,i) = 1;
            if toc(saveTimer) > save_time
                %saveData
                saveTimer = tic;
            end
            if NEGF_param.print
                disp(length(B)*(i-1)+j + "/" + length(B)*length(E));
            end
        end
    end
end
end

