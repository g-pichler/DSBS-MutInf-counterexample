%% Binary entropy function
function h=binent(x)
  h=-x.*log2(x)-(1-x).*log2(1-x);
  %% Make sure that x=1 and x=0 are handled correctly
  if (strcmp(class(x),'infsupdec'))
    I=ismember(infsupdec('0'),x) | ismember(infsupdec('1'),x);
    h(I)=union(h(I),infsupdec('0'));
  else
    h(isnan(h))=0;
  end
end
