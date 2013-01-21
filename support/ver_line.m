function ver_line(x,mrk);

%VER_LINE       ver_line([x,mrk]);
%
%               Draws a vertical dotted line along "x" on current plot
%               and leaves plot in "hold" state it was in
%
%     INPUTS:   x:      (Optional, DEF=0) The value of the horizontal axis to draw a
%                       vertical line.
%		mrk:	(Optional, DEF=':') The line type used.
%
%     OUTPUTS:          Plots line on current plot
%
%Emma     5/26/88: original
%Woodley  6/16/95: Matlab 4.1
%
%===========================================================================

if exist('x')==0,
  x = 0;                           % default to line at 0 if not given
end
if exist('mrk')==0,
  mrk = ':';                       % default to ':'
end

hold_state = get(gca,'NextPlot');  % get present hold state
YLim = get(gca,'YLim');            % get present axis limits
hold on                            % hold current plot
plot(x*ones(size(YLim)),YLim,mrk)  % draw line
hold off                           % remove hold

set(gca,'NextPlot',hold_state);    % restore original hold state
