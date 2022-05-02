function [procSample] = sample_compress(sample)
%SAMPLE_COMPRESS Summary of this function goes here
%   Detailed explanation goes here

procSample = sample;
if ~sample.compressed
    [t, method] = compress(sample.units);
    procSample.compressed = true;
    procSample.compMethod = method;
    procSample.units = t;
else
    procSample.units = decompress(sample.units, sample.compMethod);
    procSample.compressed = false;
end

