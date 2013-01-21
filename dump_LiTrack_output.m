function dump_LiTrack_output(zj,dE_Ej,fnout);

%	dump_LiTrack_output(fnout);
%
%	Function to write 2-column LiTrack output ASCII file
%	(z in meters, and dE/E in unitless fraction)
%
%	INPUTS:		zj:			The longitudinal bunch coordinate
%							of each electron [m]
%				dE_Ej:		The relative energy deviation of
%							each electron [ ]
%				fnout:		The full filename as a character
%							string

%============================================================

fid = fopen(fnout,'w');
N = length(zj);
for j = 1:N
  fprintf(fid,'%12.9e   %12.9e\n',zj(j),dE_Ej(j));
end
fclose(fid);
