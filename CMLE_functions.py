from __future__ import division
from collections import OrderedDict
import numpy as np
from adolc import condassign
from adolc import erf
from adolc import adouble 

def norm_cdf(x):
    return 0.5 + 0.5*(erf(x/np.sqrt(2.0)))






def binorm_cdf(a, b, r):
  w0 = np.array([0.04691008, 0.23076534, 0.5, 0.76923466, 0.95308992])
  w1 = np.array([0.018854042, 0.038088059, 0.0452707394, 0.038088059, 0.018854042])
  h1 = a
  h2 = b
  h12 = (h1 * h1 + h2 * h2) / 2.0
  r2 = 1.0 - r * r
  r3 = np.sqrt(r2)
  h3 = h1 * h2
  h7 = np.exp(-h3 / 2.0)
  h6 = np.abs(h1 - h2)
  h5 = h6 * h6 / 2.0
  h6 = h6 / r3
  AA = 0.5 - h3 / 8.0
  ab = 3.0 - 2.0 * AA * h5
  LH = 0.13298076 * h6 * ab * (1 - norm_cdf(h6)) - np.exp(-h5 / r2) * (ab + AA * r2) * 0.053051647
  for i in xrange(5):
    r1 = r3 * w0[i]
    rr = r1 * r1
    r2 = np.sqrt(1.0 - rr)
    LH = LH - w1[i] * np.exp(-h5 / rr) * (np.exp(-h3 / (1 + r2)) / r2 / h7 - 1.0 - AA * rr)
  
    
  h12 = condassign(h12,h1-h2,h2,h1)
  output = LH * r3 * h7 + norm_cdf(h12)
  return output







def vec2U_regr(x, regr):
    idx0 =  [np.size(regr[v],1) for v in regr.keys()]
    idx1 = np.cumsum(idx0)
    idx0 = idx1-idx0
    idx= np.vstack((idx0, idx1))
    indx = (np.arange(idx[0,0], idx[1,0]),
            np.arange(idx[0,1], idx[1,1]),
np.arange(idx[0,2], idx[1,2]),
np.arange(idx[0,3], idx[1,3]),
np.arange(idx[0,4], idx[1,4]),
np.arange(idx[0,5], idx[1,5]),
np.arange(idx[0,6], idx[1,6]))

    param = OrderedDict([('barWA', regr['barWA'].dot(x[indx[3]])),
                         ('barWB', regr['barWB'].dot(x[indx[4]])),
('bara', regr['bara'].dot(x[indx[5]])),
('VA', regr['VA'].dot(x[indx[1]])),
('VB', regr['VB'].dot(x[indx[6]])),
('SA', regr['SA'].dot(x[indx[0]])),
('CB', regr['CB'].dot(x[indx[2]])),
('sig', 1)])

    for k in param.keys():
        if not np.any(param[k]):
            param[k] = np.zeros(param[k].shape)
            
    return param

def f_jo(p, U,i):
    return norm_cdf((p*U['barWB'][i] + (1-p)*U['VB'][i] - U['CB'][i] )/(p))

def cStar(p, U,i):
    return (U['SA'][i]  - (1-p)*U['VA'][i] )/p
    
def h_jo(c, U,i):
    d1 = (U['barWA'][i] - U['bara'][i] )/(np.sqrt(2.0))
    d2 = (U['barWA'][i] - c)
    return binorm_cdf(d1, d2, 1.0/np.sqrt(2))
    
def g_jo(c, U,i):
  v1 = (c-U['barWA'][i])
  v2 = (c-U['bara'][i])
  output = (1.0 - norm_cdf(v1)*norm_cdf(v2))
  return output

  
def eqProbs(p, U):
    A = len(p)
    pC = adouble(np.zeros((A,)))
    pF = adouble(np.zeros((A,)))
    for i in xrange(A):
        c = cStar(p[i], U,i)
        g1 = g_jo(c, U,i)
        pC[i] = g1
        pF[i] = h_jo(c, U,i)/g1
    return np.hstack((p, pC, pF)).reshape((-1,3), order='F')

  
def const(p, U):
    A = len(p)
    psi = adouble(np.zeros((A,)))
    for i in xrange(A):
        c = cStar(p[i], U,i)
        g1 = g_jo(c, U,i)
        j = h_jo(c, U,i)/g1
        psi[i] = (p[i] - f_jo(j, U,i))
    return psi
    
    
  
def LL_jo(x, Y, regr):
    M = Y.shape[1]
    xP = 1.0/(1+np.exp(-x[(len(x)-M):]))
    xT = x[:(len(x)-M)]
    U = vec2U_regr(xT, regr)
    EQ = eqProbs(xP, U)
    OUT = np.hstack((1-EQ[:,1, None], EQ[:,1, None]*(1-EQ[:,0, None]), EQ[:,1, None]*EQ[:,0, None]*EQ[:,2, None], EQ[:,1, None] * EQ[:,0, None]* (1-EQ[:,2, None])))
    LL = np.sum(np.log(OUT.T) * Y)
    return -LL 
    

def const_cmle(x, Y, regr):
    M = Y.shape[1]
    xP = 1.0/(1+np.exp(-x[(len(x)-M):]))
    xT = x[:(len(x)-M)]
    U = vec2U_regr(xT, regr)
    return const(xP,U)
    