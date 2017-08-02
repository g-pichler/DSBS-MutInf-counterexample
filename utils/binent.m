%% Binary entropy function
function h=binent(x)
  h=-x.*log2(x)-(1-x).*log2(1-x);
  if ~strcmp(class(x),'infsupdec')
    h(isnan(h))=0;
  end
end
