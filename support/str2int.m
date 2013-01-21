function i = str2int(str);

%       i = str2int(str)
%
%       Converts a string to an integer.
%
%     INPUTS:   str:    The character string to convert to an integer
%
%     OUTPUTS:  i:      The integer as output

%==========================================================================

ii = find(str~=' ');
num_chars = length(ii);
if num_chars < 1
  error('STR2INT conversion error: No non-blank characters')
end

str = str(ii);

i = 0;

for j = 1:num_chars
  if ((str(j)+0)<48) | ((str(j)+0)>57),
    i = 0;
    error(['STR2INT conversion error: Non-convertible string = ' str])
  end
  i = i + (10^(num_chars-j))*(str(j)-48);
end
