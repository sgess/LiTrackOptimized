function ys = fast_smooth2(y);

[r,c] = size(y);
if r<3 | c<3
  error('fast_smooth2 does not work on matrices less than 3-by-3 - quitting.')
end
for k = 1:c
  for j = 1:r
	if k==1 & j==1
	  ys(j,k) = mean([y(j,k) y(j+1,k) y(j,k+1)]);
	elseif k==1 & j==r
	  ys(j,k) = mean([y(j,k) y(j-1,k) y(j,k+1)]);
	elseif k==c & j==1
	  ys(j,k) = mean([y(j,k) y(j+1,k) y(j,k-1)]);
	elseif k==c & j==r
	  ys(j,k) = mean([y(j,k) y(j-1,k) y(j,k-1)]);
	elseif k==1
	  ys(j,k) = mean([y(j,k) y(j-1,k) y(j+1,k) y(j,k+1)]);
	elseif k==c
	  ys(j,k) = mean([y(j,k) y(j-1,k) y(j+1,k) y(j,k-1)]);
	elseif j==1
	  ys(j,k) = mean([y(j,k) y(j,k-1) y(j,k+1) y(j+1,k)]);
	elseif j==r
	  ys(j,k) = mean([y(j,k) y(j,k-1) y(j,k+1) y(j-1,k)]);
	else
	  ys(j,k) = mean([y(j,k) y(j-1,k) y(j+1,k) y(j,k-1) y(j,k+1)]);
	end
  end
end