

function [pTime] = medFilterComplex(pTime,medFiltKernelWidth)

visPTSignal(pTime,[],[],[],1);

gpu = isa(pTime,'gpuArray');
if gpu; pTime = gather(pTime);end
pTime = medfilt1( real(pTime),medFiltKernelWidth,[],2) + ...%real channel
        1i* medfilt1(imag(pTime),medFiltKernelWidth,[],2) ; %imaginary channel
if gpu; pTime = gpuArray(pTime);end

visPTSignal(pTime,[],[],[],1);


