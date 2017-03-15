function z0=lin_interp2(x,y,z, x0,y0)
  A=[x';y';z'];
  A2=A(:,1:2)-[A(:,3),A(:,3)];
  A3=A2(1:2,:);
  if strcmp(class(A3), "infsupdec")
    singular=ismember(0, det(A3));
  else
    singular=rcond(A3) < 1e-14;
  end
  if strcmp(class(x0), "infsupdec") || strcmp(class(y0), "infsupdec") || strcmp(class(A3), "infsupdec")
    myinf=ones(size(x0))*infsupdec("0","Inf");
  else
    myinf=Inf*ones(size(x0));
  end
  if ~singular
    a=[x0';y0']-A(1:2,3);
    if length(x0) == 1
      b=A3\a;
    else
      b=inv(A3)*a;
    end
    s=A2(3,:)*b+A(3,3);
    z0=s';
  else
    z0=myinf;
  end
end
