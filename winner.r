# winner.r Copyright (c) 2018 Lloyd T. Elliott.
#
# Redistribution and use in source and binary forms, with or
# without modification, are permitted provided that the
# following conditions are met:
#
# 1. Redistributions of source code must retain the above
# copyright notice, this  list of conditions and the
# following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above
# copyright notice, this list of conditions and the following
# disclaimer in the documentation and/or other materials
# provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
# CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
# INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
# GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
# BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
# OF THE USE   OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

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
