recover.data.glmmadmb = lsmeans:::recover.data.lm

lsm.basis.glmmadmb = function (object, trms, xlev, grid)
{
  contrasts = object$contrasts
  m = model.frame(trms, grid, na.action = na.pass, xlev = xlev)
  X = model.matrix(trms, m, contrasts.arg = contrasts)
  bhat = fixef(object)
  V = vcov(object)
  misc = list()
  if (!is.null(object$family)) {
    fam = object$family
    misc$tran = object$link
    misc$inv.lbl = "response"
    if (!is.na(pmatch(fam,"binomial")))
      misc$inv.lbl = "prob"
    else if (!is.na(pmatch(fam,"poisson")))
      misc$inv.lbl = "rate"
  }
  nbasis = matrix(NA)
  dffun = function(...) NA
  list(X = X, bhat = bhat, nbasis = nbasis, V = V, dffun = dffun,
       dfargs = list(), misc = misc)
}
