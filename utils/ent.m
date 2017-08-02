%% Compute the entropy of a vector along the 1st dimension
function h=ent(x)
  t=-x.*log2(x);
  if ~strcmp(class(x), 'infsupdec')
    t(isnan(t))=0;
  end
  h=sum(t,1);
end
