from parmap import *
import numpy as np
from CMLE_sanctionsMC_support import *
from CMLE_functions import *
from rpy2.robjects import r
from rpy2.robjects.numpy2ri import numpy2ri

tempOut = r("Results <- list()")
B = 1000

i=0
hold = np.array(parmap(lambda b:MCdo(b,i), xrange(B)))
hold=numpy2ri(hold)
tempOut = r.assign("hold", hold)
tempOut = r("Results <- hold")
tempOut = r("save.image('CMLE_sanctionsMC.rdata')")

