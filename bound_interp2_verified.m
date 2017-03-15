function z0=bound_interp2_verified(x,y,z, T, x0, y0)
  i=tsearchni([x,y], T, [x0, y0]);
  assert(~isnan(i));
  t=T(i,:);
  assert(check_hull(x,y,z,t));
  z0=lin_interp2(x(t),y(t),z(t), x0,y0);
end

function i=tsearchni(X, T, A)
  if strcmp(class(X), "infsupdec")
    X1=sup(X);
  else
    X1=X;
  end
  if strcmp(class(A), "infsupdec")
    A1=sup(A);
  else
    A1=A;
  end
  
  i=tsearchn(X1, T, A1);
end

function b=check_hull(x,y,z,t)
  I=~sum( (1:length(x))' == t, 2);
  s=lin_interp2(x(t),y(t),z(t),x(I),y(I));
  if strcmp(class(s), "infsupdec")
    s=sup(s);
  end
  if strcmp(class(z), "infsupdec")
    z2=sup(z(I));
  else
    z2=z(I);
  end
  b=logical(prod(z2 <= s));
end
