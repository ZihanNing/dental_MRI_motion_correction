
function [] = visColor(co)

NColors = size(co,1);
N = 500;
X = ones([N N 3]);
mSt = makeBins(N, NColors);

for i=1:NColors
    for j=1:3
        X(mSt==i,:,j) = co(i,j);
    end
end

figure
imshow(X,[]);colormap(co)