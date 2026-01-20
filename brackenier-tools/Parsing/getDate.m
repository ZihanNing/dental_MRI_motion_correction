
function [suff, day, month, year, hour, min, sec] = getDate()

%GETDATE   gets the current date.
%   [SUFF,D,M,Y]=GETDATE()
%   ** SUFF is a string "year-month-day" to be used as a suffix (for e.g. file saving)
%   ** D is the day
%   ** M is the month
%   ** Y is the year
%
%   Yannick Brackenier 2023-01-14

%%% Get clock data
c=clock; 

%%% Split in day/month/year
day = c(3);
month = c(2);
year = c(1);
hour = c(4);
min = c(5);
sec = c(6);

%%% Create suffix
suff = sprintf('%04.f-%02.f-%02.f',year,month,day);%0 to pad with zeros rather than spaces

