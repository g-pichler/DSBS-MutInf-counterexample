classdef ibinent
  properties (Hidden, SetAccess = private)
    N;
    hs;
    xs;
    hs_int;
    xs_int;
  end
  
  methods
    function self = ibinent(N=1e5)
      pkg load interval;
      self.xs=linspace(0,1/2,N);
      self.xs_int=infsupdec(self.xs);
      self.hs=binent(self.xs);
      self.hs_int=binent(self.xs_int);
      self.N=N;
    end

    function x=binentinv(self, h)
      x=interp1(self.hs,self.xs,h);
    end

    function x=binentinv_int(self, h)

      common=subset(h, infsupdec('0','1'));
      if ~common
	h=intersect(h, infsupdec('0','1'));
      end
      
      %% start bisecting
      lo=1;
      s=self.N;
      change=true;
      while ( s > 1 || change ) && lo < self.N
	change=false;
	s=max(1,floor(s/2));
	i=lo+s;
	h1=self.hs_int(i);
	if h >= h1
	  lo=i;
	  change=true;
	end
      end
      
      hi=self.N;
      s=self.N;
      change=true;
      while ( s > 1 || change ) && hi > 1
	change=false;
	s=max(1,floor(s/2));
	i=hi-s;
	h1=self.hs_int(i);
	if h < h1
	  hi=i;
	  change=true;
	end
      end

      %% interpolate for the upper bound
      x_hi=self.lin_interp1(self.hs_int([lo,hi]), self.xs_int([lo,hi]), h);

      %% interpolation for lower bound
      x_lo=-1;
      if lo > 1
	x_lo=max(x_lo, self.lin_interp1(self.hs_int([lo-1,lo]), self.xs_int([lo-1,lo]), h));
      end
      if hi < self.N
	x_lo=max(x_lo, self.lin_interp1(self.hs_int([hi,hi+1]), self.xs_int([hi,hi+1]), h));
      end

      x=union(x_hi,x_lo);
      x=intersect(x, infsupdec('0','1/2'));
      if common
	x=infsupdec(x,'com');
      else
	x=infsupdec(x,'trv');
      end
    end
  end

  methods (Static, SetAccess = private)
    function y0=lin_interp1(X,Y,x0)
      y0=Y(1)+(x0-X(1))/(X(2)-X(1))*(Y(2)-Y(1));
    end
  end
  

end

