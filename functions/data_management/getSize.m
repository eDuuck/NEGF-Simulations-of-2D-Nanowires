function totSize = getSize(obj,depth)
%GETSIZE Summary of this function goes here
%   Detailed explanation goes here
if ~exist("depth","var")
    depth = inf;
end
totSize = 0;
if isstruct(obj)
    props = fieldnames(obj);
elseif iscell(obj)
    for i = 1:length(obj)
        currentProperty = obj{i};
        totSize = totSize + getSize(currentProperty,depth);
    end
    return
else
    props = properties(obj);
end
if isempty(props) || depth <= 0
    s = whos("obj");
    totSize = s.bytes;
    return
else
    for i = 1:length(props)
        currentProperty = obj.(char(props(i)));
        totSize = totSize + getSize(currentProperty,depth-1);
    end
end

