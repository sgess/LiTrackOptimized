function output_string = lowcase(input_string)

%  Function result = lowcase(input) 
%
% Converts the upper case letters into lower case.

% find the size so we can loop through even a douple-dimension 
% string

[axx,bxx]=size(input_string);

% allocate a string of the same dimensions as the input

output_string = input_string;

% loop through all characters, testing for case.  Note that a 
% douple dimension array can be treated as a single dimension.

for i = 1:axx*bxx
  if output_string(i) >= 'A' & output_string(i) <= 'Z',
    output_string(i) = output_string(i) + 'a' - 'A';
  end
end

% set the returned string back into string syntax

output_string = setstr(output_string);
