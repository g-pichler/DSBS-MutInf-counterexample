function q=hyperplane3(x,y,z)
  A=[x, y, z]';
  a=A(:,1);
  b=A(:,2);
  c=A(:,3);
  d=cross(interv_mean(a-c),interv_mean(b-c));

  if d(3) > 0
    d=-d;
  end
  
  k=d'*[a,b,c];
  k_inf=inf(min(k));

  q=[d; k_inf];
end

function b=interv_mean(a)
  b=(sup(a)+inf(a))/2;
end
