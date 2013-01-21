function hor_line(y,mrk);

%HOR_LINE       hor_line([y,mrk]);
%
%               Draws a horizontal dotted line along "y" on current plot
%               and leaves plot in "hold" state it was in
%
%     INPUTS:   y:      (Optional, DEF=0) The value of the vertical axis to draw a
%                       horizontal line.
%		mrk:	(Optional, DEF=':') The line type used.
%     OUTPUTS:          Plot line on current plot
%
%Emma     5/26/88: original
%Woodley  6/16/95: Matlab 4.1
%
%===========================================================================

if exist('y')==0,
  y = 0;                           % default to line at 0 if not given
end
if exist('mrk')==0,
  mrk = ':';                       % default to ':'
end

hold_state = get(gca,'NextPlot');  % get present hold state
XLim = get(gca,'XLim');            % get present axis limits

hold on                            % hold current plot
plot(XLim,y*ones(size(XLim)),mrk)  % draw line
hold off                           % remove hold

set(gca,'NextPlot',hold_state);    % restore original hold state
