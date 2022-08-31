function output = NEGF_map(sample,E,B,result,saveTime,errorMarg,rate,it_lim,reduce)
%MC_NEGF Summary of this function goes here
%   Detailed explanation goes here
if ~exist("errorMarg","var")
    errorMarg = 1e-6 * min(sample.getUnits,[],'all');
end
if ~exist("rate","var")
    rate = 0.5;
end
if ~exist("it_lim","var")
    it_lim = 50;
end
if ~exist("reduce","var")
    reduce = false;
end
if ~exist("saveTime","var")
    saveTime = Inf;
end
newRes = true;
if exist("result","var")
    if isstruct(result)
        output = result;
        newRes = false;
    end
end
if newRes
    output = struct("E",E,...
    "B",B,"completed", zeros(length(B),length(E)));
    output.NEGF_results = cell(length(B),length(E));
end
tic
disp(0 + "/" + length(B)*length(E));
for i = 1:length(E)
    for j = 1:length(B)
        if ~output.completed(j,i)
            r0 = 0;
            if j ~= 1
                r0 = output.NEGF_results{j-1,i};
            end
            output.NEGF_results{j,i} = NEGF(sample,E(i),B(j),errorMarg,...
                            rate,it_lim,reduce,r0);
            output.completed(j,i) = 1;
            if toc > saveTime
                %saveData
                tic
            end
            disp(length(B)*(i-1)+j + "/" + length(B)*length(E));
        end
    end
end
end

