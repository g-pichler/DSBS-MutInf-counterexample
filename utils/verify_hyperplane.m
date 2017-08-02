function b=verify_hyperplane(q,x,y,z)
  a=[x, y, z]*q(1:3);
  b=all(precedes(q(4), a));
end
