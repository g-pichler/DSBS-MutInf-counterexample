function z=z_interp2(q,x,y)
  z = ( q(4) - (q(1:2)'*[x;y]) ) / q(3);
end
  
