        function s = prompt(str,options,deflt_option)

%       s = prompt(str,options[,deflt_option]);
%
%	Prompts user for input of ONE character from "options" string.  
%	It uses only the first character of the user input, converted to 
%	lower case, and returns this character which will then match 
%	at least one of the characters in the "options" string.
%	When in batch mode, it automatically chooses the one 
%	character string "batch_choice" without prompting.  Each 
%	character in the "options" string is an acceptable
%	choice.  If the user types a character not in the 
%	options string, the function displays a message and prompts 
%	again (in interactive mode).
%
%     INPUTS:   str:            The prompt string displayed to the screen
%                               as the question to be answered (a list of
%                               options and a colon are tacked on at the
%                               end).
%
%               options:        A character string of options (up to 5 
%                               characters) which are possible choices 
%                               the user has.
%
%		deflt_option:	(Optional,DEF=no default) A character which
%				is the default selection when a <CR> is hit
%				(will use this choice in batch mode).
%
%   OUTPUTS:    s:              The character finally chosen by the user 
%                               (or by batch auto) which has been converted
%                               to lower case, is guaranteed to be only one
%                               character in length, and matches at least 
%                               one of the characters in the "options"
%                               string.
%
%   EXAMPLE:	>>opts = ['yn'];
%               >>s = prompt('Plot this again',opts);
%
%                 Plot this again [y/n]: 
%
%               (Now the user may only enter 'y', or 'Y', or 'n', or 'N')

%===============================================================================

n = length(options);

if n <= 1
  error('less than 2 options makes no sense')
end

if n > 5
  error('PROMPT works with up to 5 options only - see "menus"')
end

options = lowcase(options);

%if ~isitaterm('sys$output:')
%  if ~exist('deflt_option')
%    deflt_option = options(1);
%  else
%    deflt_option = lowcase(deflt_option(1));
%    ii = find(deflt_option == options);
%    if length(ii) == 0
%      error('"default_option" must be one of the provided "options"')
%    end
%  end
%  s = lowercase(deflt_option(1));
%  return
%else
  if ~exist('deflt_option')
    no_default = 1;
  else
    deflt_option = lowcase(deflt_option(1));
    ii = find(deflt_option == options);
    if length(ii) == 0
      error('"default_option" must be one of the provided "options"')
    end
    no_default = 0;
  end
%end

opts = ' (';
for j = 1:(n-1)
  opts = [opts options(j) '/'];
end
opts = [opts options(n) ')'];
if no_default
  opts = [opts ': '];
else
  opts = [opts ' [' deflt_option ']: '];
end
str = [str opts];

while 1
  s = input(str,'s');
  ns = length(s);
  if (ns==0) & (~no_default)
    s = deflt_option;
    return
  end
  if (ns==0) & (no_default)
    disp2('No default answer... select from the following:')
    disp(options(:))
  end  
  if ns ~= 0
    s = lowcase(s(1));
    for j = 1:n
      if s == options(j)
        return
      end
    end
    disp2('Select only from the following:')
    disp(options(:))
  end
end
