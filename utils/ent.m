%% Compute the entropy of a vector along the 1st dimension
function h=ent(x)
  t=-x.*log2(x);
  %% Make sure that x=1 and x=0 are handled correctly
  if (strcmp(class(x),'infsupdec'))
    I=ismember(infsupdec('0'),x) | ismember(infsupdec('1'),x);
    t(I)=union(t(I),infsupdec('0'));
  else
    t(isnan(t))=0;
  end
  h=sum(t,1);
end
