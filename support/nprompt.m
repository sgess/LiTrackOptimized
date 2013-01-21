        function num = nprompt(str,deflt_num,min_num,max_num)

%       num = nprompt(str[,deflt_num,min_num,max_num]);
%
%	Prompts user for input of a number (integer or real).  A
%       default choice may be selected by passing a default number
%       "deflt_num" to the function which is selected when the user
%       simply hits return.  In batch mode, it automatically chooses
%       "deflt_num" without a prompt.
%
%     INPUTS:   str:            The prompt string displayed to the screen
%                               as the question to be answered (the
%                               default number and a colon is appended to
%                               the end).
%
%		deflt_num:	(Optional,DEF=no default) A number which
%				is the default selection when a <CR> is hit
%				(will use this choice in batch mode).
%
%               min_num:        (Optional,DEF=+inf) The minimum acceptable
%                               number to enter.
%
%               max_num:        (Optional,DEF=-inf) The maximum acceptable
%                               number to enter.
%
%   OUTPUTS:    num:            The number finally chosen by the user 
%                               (or by batch auto) where 
%                                       min_num <= num <= max_num
%                               if "max_num" and "min_num" are provided.
%
%   EXAMPLE:	>>num = nprompt('Enter frequency in MHz',8.5);
%
%                 Enter frequency in MHz [8.5]: 
%
%               (Now the user may enter any number, or <CR> = 8.5)

%===============================================================================

if nargin == 3
  error('Need both "max_num" and "min_num" arguments')
end

bad = 1;

if ~exist('deflt_num')
  str1 = [str ': '];
  while bad
    num = input(str1);
    if length(num) ~= 0
      bad = 0;
    end
  end
else
  if exist('max_num')
    if deflt_num > max_num
      error('default selection > maximum allowable selection')
    elseif deflt_num < min_num
      error('default selection < minimum allowable selection')
    end
    str1 = [str sprintf(' (%g to %g) [%g]: ',min_num,max_num,deflt_num)];
    while bad
      num = input(str1);
      if length(num) == 0
        num = deflt_num;
      end
      if num > max_num
        disp(sprintf('Entry out of range (>%g)',max_num))
      elseif num < min_num
        disp(sprintf('Entry out of range (<%g)',min_num))
      else
        bad = 0;
      end
    end
  else
    str1 = [str sprintf(' [%g]: ',deflt_num)];
    num = input(str1);
    if length(num) == 0
      num = deflt_num;
    end
  end
end
