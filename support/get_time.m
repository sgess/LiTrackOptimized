function date_time = get_time

%
%       Function to return present date/time in character ASCII
%       string format (e.g. 03-APR-1989 17:42:15).  No input arguments.

%============================================================================


months = ['JAN'; 'FEB'; 'MAR'; 'APR'; 'MAY'; 'JUN'; ...
         'JUL'; 'AUG'; 'SEP'; 'OCT'; 'NOV'; 'DEC'];

t = fix(clock);

day = int2str(t(3));
if length(day) == 1
  day = ['0' day];
end

month = months(t(2),:);

year  = int2str(t(1));

hour  = int2str(t(4));
if length(hour) == 1
  hour = ['0' hour];
end

min   = int2str(t(5));
if length(min) == 1
  min = ['0' min];
end

sec   = int2str(t(6));
if length(sec) == 1
  sec = ['0' sec];
end

date_time = [day '-' month '-' year ' ' hour ':' min ':' sec];
