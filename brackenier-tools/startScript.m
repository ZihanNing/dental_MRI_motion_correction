
function [currentDir] = startScript()

clc; close all; clear all; 
currentDir = mfilename('fullpath');
cd(fileparts(currentDir));