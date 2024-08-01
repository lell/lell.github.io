library(neldermead)

# Winner takes beta, se, N and alpha
# Winner returns a corrected beta (Xiao & Boehnke, 2011)

winner = function(beta, se, N, alpha = 0.05, tol = 1e-8) {
  
  # We'll consider positive beta, and add the sign back at the end

  sig = sign(beta)
  beta = beta * sig
  T = beta / se
  se2 = se * se

  # For sub-threshold T, the likelihood is 0, and so beta1 = 0 maximises it

  thresh = pt(alpha / 2.0, N - 2) 
  if (T < thresh) {
    return(0)
  }

  # Curry the objective function

  obj = function(beta1) {
    numer = dnorm(beta, beta1, se2)
    denom = 1 - pt(alpha / 2.0, N - 2, ncp = beta1 / se)
    result = numer / denom

    # We'll minimise the objective function
    # So, flip the sign, so we maximise likelihood

    result = -result
    return(result)
  }

  # Minimise objective function using Nelder-Mead, initialise empirically

  result = fminsearch(fun = obj, x0 = beta)
  result = neldermead.get(result, "xopt")[1]

  # add sign to result

  result = sig * result
  return(result)
}
