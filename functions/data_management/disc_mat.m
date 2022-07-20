function result = disc_mat(A,method,byteSize)
%DISCRETIZE Summary of this function goes here
%   Detailed explanation goes here
A = full(A);

if ~exist('byteSize','var')
    byteSize = 1;
end
if ~exist('method','var')
    method = "linear";
end

if lower(method) == "linear"
    realVal = [min(real(A),[],'all'),max(real(A),[],'all')];
    imagVal = [min(imag(A),[],'all'),max(imag(A),[],'all')];
    range = [realVal;imagVal];
    switch byteSize
        case 1
            B = uint8((real(A)-realVal(1))*255 /diff(realVal) + ...
                (imag(A)-imagVal(1))*255i/diff(imagVal));
        case 2
            B = uint16((real(A)-realVal(1))*65535 /diff(realVal)+ ...
                (imag(A)-imagVal(1))*65535i/diff(imagVal));
        case 4
            B = uint32((real(A)-realVal(1))*(2^32-1)/diff(realVal) + ...
                1i*(imag(A)-imagVal(1))*(2^32-1)/diff(imagVal));
        case 8
            B = uint64((real(A)-realVal(1))*(2^64-1)/diff(realVal) + ...
                1i*(imag(A)-imagVal(1))*(2^64-1)/diff(imagVal));
        case deafult
            error('byteSize needs to be of size{1,2,4,8}')
    end
    result = struct('matrix',B,'range',range,'byteSize',byteSize,'method',"Linear");
elseif lower(method) == "non-linear"
    [realRange,realVals] = non_lin_dest(real(A),byteSize);
    [imagRange,imagVals] = non_lin_dest(imag(A),byteSize);
    matrix = uint8(realVals + 1i* imagVals);
    result = struct('matrix',matrix,'range',[realRange;imagRange],...
        'byteSize',byteSize,'method',"non-linear");
else
    error("Specified method not supported.")
end
end

function [range,values] = non_lin_dest(data,byteSize)
range = zeros(1,2^(byteSize*8));
data1D = sort(reshape(data,[],1));

range(1) = data1D(1);
range(end) = data1D(end); 
amplitude = max(abs(range(1)),range(end));


[destribution,vals] = histcounts(data/amplitude,2^(byteSize*8+6));
int = cumsum(destribution);
CRM = 30;
p1 = find(int>int(end)/CRM,1,'first');
p2 = find(int<29*int(end)/CRM,1,'last');

EPR = 100;
EPM = EPR;

lowPoints = 1;
k = 2;
while lowPoints < 2^(byteSize*8-3) && int(k) < int(end)/CRM
    if data1D(k)-data1D(lowPoints) > amplitude / EPM || ...
            data1D(k)-data1D(k-1) > amplitude / EPR
        lowPoints = lowPoints + 1;
        range(lowPoints) = data1D(k);
    end
    k = k+1;
end
lp = data1D(k+1);
highPoints = 1;
k = 1;
while highPoints < 2^(byteSize*8-3) && int(k) < (CRM-1)*int(end)/CRM
    if range(end-highPoints+1)-data1D(end-k) > amplitude / EPM || ...
        data1D(end-(k-1))-data1D(end-k) > amplitude / EPR 

        highPoints = highPoints + 1;
        range(end-highPoints+1) = data1D(end-k);
    end
    k = k+1;
end
hp = data1D(end-k-2);
if ~isempty(find(data1D==0,1))
    closerange = sort([vals(floor(linspace(p1,p2,2^(8*byteSize-1)-1))),0]);
else
    closerange = vals(floor(linspace(p1,p2,2^(8*byteSize-1))));
end
lowrange = linspace(lp,vals(p1),2^(byteSize*8-2)-lowPoints+2);
highrange = linspace(vals(p2),hp,2^(byteSize*8-2)-highPoints+2); 
range = [range(1:lowPoints), lowrange(2:end-1), closerange, highrange(2:end-1),range(end-highPoints+1:end)];
values = discretize(data,range);
end