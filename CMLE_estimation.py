import shutil
import os
from time import time
from CMLE_functions import *
import adolc 
import scipy.optimize as opt
import pyipopt    
from rpy2.robjects import r
from rpy2.robjects.numpy2ri import numpy2ri
import numpy as np
from collections import OrderedDict
import pickle as pk
from scipy.stats import norm
import warnings

user = os.path.expanduser("~")

r('source("CMLE_estimation_support.R")')
regrMat = np.array((r("do.call(cbind, regr)")))
regr = OrderedDict([('SA',regrMat[:,:4]),('VA',regrMat[:,4:6]),('CB',regrMat[:,6:11]),('barWA',regrMat[:,11:15]),('barWB',regrMat[:,15:18]),('bara',regrMat[:,18:20]),('VB',np.ones((regrMat.shape[0],0)))])
Y = np.array((r("Y")))
x0 = np.array((r("unname(x0)")))
xL = np.array((r("unname(xL)")))
   
def fll(x):
    return LL_jo(x, Y, regr)

def heq(x):
    return const_cmle(x, Y, regr)

ccd = os.getcwd()
if not os.path.exists(user + '/Documents/adolcSanctions'):
    os.makedirs(user + '/Documents/adolcSanctions')

os.chdir(user + '/Documents/adolcSanctions')



adolc.trace_on(1)
ax = adolc.adouble(np.zeros(len(x0)))
adolc.independent(ax)
ay = fll(ax)
adolc.dependent(ay)
adolc.trace_off()

# trace constraint function
adolc.trace_on(2)
ax = adolc.adouble(x0)
adolc.independent(ax)
ay = heq(ax)
adolc.dependent(ay)
adolc.trace_off()

npar = len(x0)

def adFun(x):
    return adolc.function(1, x)

def grFun(x): 
    return adolc.gradient(1, x)

def  const_adolc(x):
    return adolc.function(2,x)

def jac_adolc(x):
    return adolc.jacobian(2,x)


def lagrangian(x, lagrange, obj_factor):
    return  obj_factor*fll(x) + np.dot(lagrange, heq(x))



####Jacobian####


## initalize it
class jac_c_adolc:
    
    def __init__(self, x):
        #options = np.array([1,1,0,0],dtype=int)
        options = None
        result = adolc.colpack.sparse_jac_no_repeat(2,x,options)
        
        self.nnz  = result[0]     
        self.rind = np.asarray(result[1],dtype=int)
        self.cind = np.asarray(result[2],dtype=int)
        self.values = np.asarray(result[3],dtype=float)
        
    def __call__(self, x, flag, user_data=None):
        if flag:
            return (self.rind, self.cind)
        else:
            result = adolc.colpack.sparse_jac_repeat(2, x, self.nnz, self.rind,
                self.cind, self.values)
            return result[3]

##### create the function
Jac_c_adolc = jac_c_adolc(x0)


###Hessian
    
# trace lagrangian function
adolc.trace_on(3)
ax = adolc.adouble(x0)
adolc.independent(ax)
ay = lagrangian(ax, xL, np.array([1.0]))
adolc.dependent(ay)
adolc.trace_off()

M = Y.shape[1]
nreal = npar-M
given = {'rind': np.concatenate((np.kron(np.arange(nreal), np.ones(npar,dtype='int')), np.arange(nreal, npar))),
         'cind': np.concatenate((np.tile(np.arange(npar), nreal), np.arange(nreal, npar)))}
mask = np.where(given['rind'] <= given['cind'])
given['rind'] = given['rind'][mask]
given['cind'] = given['cind'][mask]


def hessLag_adolc(x, lagrange, obj_factor, flag, user_data=None):
    if flag:
        result = (given['rind'], given['cind'])
    else:
         result = np.ravel(adolc.hessian(3, x)[given['rind'],given['cind']], order="C")
    return result


H2 = hessLag_adolc(x0, xL, 1.0, False)
H2a = hessLag_adolc(x0, xL, 1.0, True)
nnzh = len(given['rind'])

##Optimization
#PRELIMS: other things to pass to IPOPT
nvar = len(x0) #number of variables in the problem
x_L = np.array([-np.inf]*nvar, dtype=float) #box contraints on variables (none)
x_U = np.array([np.inf]*nvar, dtype=float)
 
#PRELIMS:define the (in)equality constraints
ncon = heq(ax).shape[0] #number of constraints
g_L = np.array([0]*ncon, dtype=float) #constraints are to equal 0
g_U = np.array([0]*ncon, dtype=float) #constraints are to equal 0


#PRELIMS: define the number of nonzeros in the jacobian 
val = Jac_c_adolc(x0, False) 
nnzj = len(val)            

  
# create the nonlinear programming model
nlp2 = pyipopt.create(
nvar, 
x_L,
x_U,
ncon,
g_L,
g_U,
nnzj,
nnzh,
adFun,
grFun,
const_adolc,
Jac_c_adolc,
hessLag_adolc
)

nlp2.num_option('expect_infeasible_problem_ctol', 1e-15)
nlp2.int_option('max_iter', 5000)

nlp2.num_option('dual_inf_tol', 1e-5)
nlp2.num_option('constr_viol_tol', 1e-5)
nlp2.num_option('tol', 1e-6)


out = -np.ones(6)
xin = x0
count =0
while out[5] == -1: # iter limit exceed
    out = nlp2.solve(xin)
    xin = out[0]
    count+=1


# free the model
nlp2.close()


os.chdir(ccd)
shutil.rmtree(user + '/Documents/adolcSanctions')
output = out[0]
lam = out[3]
code = out[5]
val = out[4]

print "Estimates"
print np.round(output[:20],2)

output=numpy2ri(output)
value=numpy2ri(np.array([val]))

r.assign("output", output)
r.assign("value", value)
r("save(list=c('output',  'value'), file='CMLE_estimation_output.rdata')")

warnings.warn("End of file. Press enter if the system hangs here.")