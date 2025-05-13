def binomial(n):
  '''
  It computes binomial coefficients.

  Parameters:
  n (int): it's a nonnegative integer.
  
  Returns:
  res (list): res[i][j] is the binomial coefficient C(i,j), with 0<=j<=i<=n.
  '''
  
  res = [[1]]
  for i in range(1,n+1):
    si1s = res[-1]
    sis = [1]
    for k in range(1,i):
      sis.append(si1s[k]+si1s[k-1])
    sis.append(1)
    res.append(sis)
  return res

def stirling(n):
  '''
  It computes Stirling numbers of the second kind.

  Parameters:
  n (int): it's a nonnegative integer.
  
  Returns:
  res (list): res[i][j] is the Stirling number of the second kind S(i,j), with 0<=j<=i<=n.
  '''

  res = [[1]]
  for i in range(1,n+1):
    si1s = res[-1]
    sis = [0]
    for k in range(1,i):
      sis.append(k*si1s[k]+si1s[k-1])
    sis.append(1)
    res.append(sis)
  return res

nmax = 4              # It indicates until which moment one wishes to compute
lnu = np.log(1.01)    # It indicates the value of log(u)
lnd = np.log(1/1.01)  # It indicates the value of log(d)
t = 100
N = 1000
Nu = 500
Nd = 300

# Changing the previous parameters, the following script prints a list
# with the moments of the logarithmic return, computed using the
# deduced formula: moments[i] is the i-th moment

s = stirling(nmax)
b = binomial(nmax)
js = []
for j in range(nmax+1):
  js.append((-1)**j*(1-j/N)**t)
lnus = [1]
lnds = [1]
Nus = [1]
Nds = [1]
for j in range(nmax):
  lnus.append(lnus[-1]*lnu)
  lnds.append(lnds[-1]*lnd)
  Nus.append(Nus[-1]*(Nu-j))
  Nds.append(Nds[-1]*(Nd-j))

moments = [1]
for n in range(1,nmax+1):
  s0 = 0
  for nu in range(n+1):
    c1 = b[n][nu]*lnus[nu]*lnds[n-nu]
    s1 = 0
    nd = n-nu
    for ku in range(nu+1):
      c2 = Nus[ku]*s[nu][ku]
      s2 = 0
      for kd in range(nd+1):
        c3 = Nds[kd]*s[nd][kd]
        s3 = 0
        k = ku+kd
        for j in range(k+1):
          s3 += js[j]*b[k][j]
        s2 += c3*s3
      s1 += c2*s2
    s0 += c1*s1
  moments.append(s0)

print(moments)
