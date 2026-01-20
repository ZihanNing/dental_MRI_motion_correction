

function [memAvailable] = getGPUMem()

h = gpuDevice;

memAvailable = h.AvailableMemory;
memTotal = h.TotalMemory;

fprintf('Available GPU memory: %.2f GB / %.2f GB\n',memAvailable/1e9,memTotal/1e9);