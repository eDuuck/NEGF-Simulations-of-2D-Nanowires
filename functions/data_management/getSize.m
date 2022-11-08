function totSize = getSize(obj,depth)
%GETSIZE returns the data size of an object. If the object is a structure
%and contains pointers to other objects, these will also be included.
%   GETSIZE(obj) returns the byte size of obj. If obj is a structure or
%   cell array then every element or property of the obj will be also be
%   included. If one of the elemnts or properties are a pointer, then
%   GETSIZE will continue to look into that pointed data object. The
%   algorithm breaks automatically when reaching a depth of 5 or if the
%   object isn't a structure/cell array.
%
%   GETSIZE(obj,depth) will perform the same action as GETSIZE(obj) but the
%   depth of the search algorithm breakpoint can be specified. If a
%   structure points to another structure then this other structure can be
%   ignored by choosing an appropriete depth. In the case of NEGF_result,
%   where the result contains a pointer to the sample that has been
%   simulated, the sample can be ignored in this algorithm by calling
%   GETSIZE(NEGF_result, 0).

if ~exist("depth","var")
    depth = 5;
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

