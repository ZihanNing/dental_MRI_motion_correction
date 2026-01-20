
function [RFFlag, RFType] = convertTWIX_ucExcitMode(mode)

RFFlag = 1;

if mode==1
    RFType='Slice-Selective';
    
elseif mode==2
    RFType='Volume-Selective';
    
elseif mode==4
    RFType='Adiabatic';
    
elseif mode==8
    RFType='Multi-Slab-Selective';
    
elseif mode==10
    RFType='Slice-Selective-PE';
    
elseif mode==20
    RFType='Pulse-Standard';
    
elseif mode==40
    RFType='Zoomed';
    
elseif mode==80
    RFType='User-Defined';
    
elseif mode==100
    RFType='Tone';
    
else
    RFType='unknown';
end