%% Binary convolution
function c=star(a,b)
  c=a.*(1-b)+(1-a).*b;
end
