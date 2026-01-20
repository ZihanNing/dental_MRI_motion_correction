
function [InversionFlag, InversionType] = convertTWIX_ucInversion(mode)

InversionFlag = mode~=4;
if InversionFlag==0; InversionType='NA';return;end

if mode==1
    InversionType='Slice-Selective';
    
elseif mode==2
    InversionType='Volume-Selective';
    
elseif mode==8
    InversionType='T2-Selective';
    
elseif mode==10
    InversionType='Volume-Selective-Dir';
    
elseif mode==20
    InversionType='Slice-Selective-Dir';
    
elseif mode==40
    InversionType='Volume-Selective-Tuned';
    
elseif mode==80
    InversionType='Volume-Selective-T2-Prepared';
    
elseif mode==100
    InversionType='Volume-Selective-HighBW';
    
else
    InversionType='unknown';
end