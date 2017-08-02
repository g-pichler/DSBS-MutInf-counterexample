function L=upper_concave_envelope3_tri(x,y,z)
  xmin=min(x);
  xmax=max(x);
  ymin=min(y);
  ymax=max(y);
  z0=min(z)-1;
  delta=1e-14;

  x1=[xmin; xmax; xmin; xmax; x];
  y1=[ymin; ymin; ymax; ymax; y];
  z1=[z0;z0;z0;z0;z];
  X=[x1,y1,z1];
  
  K=convhulln(X);

  M=min(K,[],2);
  K=K(M>4,:);
  K=K-4;
  lK=size(K)(1);

  %% Remove degenerate triangles
  b=logical(ones(lK,1));
  for i=1:lK
    t=K(i,:);
    if sum(abs(x(t)-xmin)) <= delta || sum(abs(x(t)-xmax)) <= delta ...
       || sum(abs(y(t)-ymin)) <= delta || sum(abs(y(t)-ymax)) <= delta
      b(i)=0;
    end
  end
  L=K(b,:);
end
