
function [rec] = blockPreProcessing(rec)

rec.Alg.parXT.PT.signalUsage.eigTh = onesL(rec.Alg.parXT.PT.signalUsage.eigTh);
rec.Alg.parXT.PT.signalUsage.orderPreProcessing = [];
rec.Alg.parXT.PT.signalUsage.referencePhase = 0;
rec.Alg.parXT.PT.PTFilter.medFiltKernelWidthms =0;
rec.Alg.parXT.PT.PTFilter.golayFiltKernelWidthms = 0;