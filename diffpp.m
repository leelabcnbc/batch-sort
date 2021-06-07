function dpp = diffpp(pp)
% DIFFPP differentiate a pp function%
%dpp = diffpp(pp)
%
% returns the first derivative of the spline in  PP .
% 2 oct 95 cb
[breaks,coefs] = unmkpp(pp);
[l,k] = size(coefs);
if k==1,
   dpp = mkpp(breaks,zeros(1,l));
else   
   [k-1:-1:1];
   ans(ones(l,1),:).*coefs(:,1:k-1);
   dpp = mkpp(breaks,ans(:));
end