
clc
clear all
close all

t = linspace(0,2*pi,50);
fun = @(x) func(x,t);

rng(22)
xTrue = plugNoise(ones([1 2]),1,0)+[.1 pi];
yTrue = fun(xTrue);
y = yTrue + .0005*plugNoise(yTrue,0);
%y=yTrue;
x0 = xTrue + .01*plugNoise(xTrue,1);
%x0 = [];

tic
yFit = zerosL(y);
for i=1:size(y,1); yFit = dynInd(yFit,i,1,curveFit(dynInd(y,i,1), fun, x0(1,:), [], [], [], [], 2) );end
%[yFit, x] = curveFit(y, fun, x0, [], [], [], [], 2);
toc

for i=1:size(yFit,1)
    figure(100)
    clf;
    plot(t,abs(yTrue(i,:)));hold on
    plot(t,abs(yFit(i,:)));hold on
    plot(t,abs(y(i,:)));
    axis([-inf inf -inf inf])
    pause(.1)
end

function [f] = func(x,t)
 %f=x(1,1).*exp(-(t-x(1,2)).^2./(2*x(:,3).^2) ); 
 f = 1i./(pi*x(:,1)) .*(x(:,1).^2./((t-x(:,2)).^2+x(:,1).^2));%Lorentzian
end