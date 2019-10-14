import shutil
import os
from time import time
from CMLE_functions import *
import adolc 
import scipy.optimize as opt
import pyipopt    
import rpy2.robjects as robj #for loading the .rdata file
import numpy as np
from collections import OrderedDict


def MCdo(b, i):
    user = os.path.expanduser("~")

    tempOut = robj.r("rm(list=ls())")
    tempOut = robj.r("i <- %i"%(i+1))
    tempOut = robj.r("b <- %i"%(i*(1000)+(b+1) ))
    tempOut = robj.r("source('CMLE_unstable_support.r')")
    regrMat = np.array((robj.r("do.call(cbind, regr)")))
    regr = OrderedDict([('SA',np.ones((regrMat.shape[0],0))),('VA',regrMat[:,0:1]),('CB',np.ones((regrMat.shape[0],0))),('barWA',regrMat[:,1:2]),('barWB',regrMat[:,2:4]),('bara',regrMat[:,4:5]),('VB',regrMat[:,5:6])])
    Y = np.array((robj.r("Y")))
    x0 = np.array((robj.r("unname(x0)")))
    xL = np.array((robj.r("unname(xL)")))
    tPL = np.array((robj.r("unname(out.2step$time )")))
    
    def fll(x):
        return LL_jo(x, Y, regr)
    
    def heq(x):
        return  const_cmle(x, Y, regr)
    np.random.seed(b)
    while np.isinf(fll(adouble(x0)).val):
        x0 = np.random.uniform(size=len(x0))
    
    ccd = os.getcwd()
    if not os.path.exists(user + '/Documents/adolc%i_%i'%(i, b)):
        os.makedirs(user + '/Documents/adolc%i_%i'%(i, b))
    
    os.chdir(user + '/Documents/adolc%i_%i'%(i, b))
    
    
    
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

    #Jacobian
    
    
    #### initalize it
    class jac_c_adolc:
        
        def __init__(self, x):
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
    nlp2.int_option('max_iter', 100)
    nlp2.num_option('dual_inf_tol', 1e-3)
    nlp2.num_option('constr_viol_tol', 1e-3)
    nlp2.num_option('tol', 1e-6)
    nlp2.int_option('print_level', 0)
    
    t1 = time()
    out = nlp2.solve(x0)
    t1 = time()-t1
    t1 = t1 + tPL
    
    # free the model
    nlp2.close()

    
    os.chdir(ccd)
    shutil.rmtree(user + '/Documents/adolc%i_%i'%(i, b)) 
    output = np.concatenate((out[0][:6], t1, np.array([out[5]])))
    return output
