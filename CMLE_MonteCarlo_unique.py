from parmap import *
import numpy as np
from CMLE_MonteCarlo_unique_support import *
from CMLE_functions import *
from rpy2.robjects import r
from rpy2.robjects.numpy2ri import numpy2ri
import warnings

r("Results <- list()")
B = 1000


for i in xrange(20):     
    hold = np.array(parmap(lambda b:MCdo(b,i), xrange(B)))
    hold=numpy2ri(hold)
    r.assign("hold", hold)
    tempOut = r("Results[[%i]] <- hold"%(i+1))
    tempOut = r("save.image('CMLE_unique.rdata')")

warnings.warn("End of file. Press enter if the system hangs here.")
