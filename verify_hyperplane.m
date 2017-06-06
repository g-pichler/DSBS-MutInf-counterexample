function b=verify_hyperplane(q,x,y,z)
  a=q(1:3)'*[x';y';z'];
  b=all(precedes(q(4), a));
end
