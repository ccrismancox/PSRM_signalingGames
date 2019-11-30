from parmap import *
import numpy as np
from CMLE_unstable_support import *
from CMLE_functions import *
from rpy2.robjects import r
from rpy2.robjects.numpy2ri import numpy2ri

tempOut = r("Results <- list()")
B = 1000


for i in xrange(11):     
    hold = np.array(parmap(lambda b:MCdo(b,i), xrange(B)))
    hold=numpy2ri(hold)
    r.assign("hold", hold)
    tempOut = r("Results[[%i]] <- hold"%(i+1))
    tempOut = r("save.image('CMLE_unstable.rdata')")

