function y=dist(w,p)

% DIST Distances between vectors.
%        DIST(W,P)
%          W - SxR matrix of rows vectors.
%          P - RxQ matrix of column vectors.
%        Returns SxQ matrix of vectors distances.
%
%        EXAMPLE: w = [1 2 1];
%                p = [1; 2.1; 0.9];
%                 n = dist(w,p)       

[s,r] = size(w);
[r2,q] = size(p);

if (r ~= r2), error('Matrix sizes do not match.'),end

y = zeros(s,q);

if r == 1
  for i=1:s
    x = w(i,:)'*ones(1,q);
    y(i,:) = abs(x-p);
  end
else
  for i=1:s
    x = w(i,:)'*ones(1,q);
    y(i,:) = sum((x-p).^ 2) .^ 0.5;
  end
end           