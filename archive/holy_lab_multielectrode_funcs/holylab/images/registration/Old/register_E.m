function err = register_E(psiF,psig,w,sqrtdetJ,spacing,lambda)
  isnz = w(:) > 0;
  err = sum(w(isnz) .* (psiF(isnz) - psig(isnz)).^2);
  fprintf('data error: %g, ',err);
  if any(lambda)
    %keepIndex = ~isnan(sqrtdetJ);
    %err = err + lambda(1) * sum((sqrtdetJ(keepIndex(:)) - 1).^2);
    %err = err + lambda(1) * nanmean((sqrtdetJ(:) - 1).^2);
    err = err + lambda(1) * mean((sqrtdetJ(:) - 1).^2);
    if (lambda(2) ~= 0)
      sz = size(psiF);
      dimKeep = sz > 1;
      for i = 1:length(sz)
        if (dimKeep(i))
          dJ = diff(sqrtdetJ,1,i);
          %err = err + lambda(2) * nanmean(dJ(:).^2)/spacing(i)^2;
          err = err + lambda(2) * mean(dJ(:).^2)/spacing(i)^2;
        end
      end
    end
  end
  fprintf('total error %g\n',err);