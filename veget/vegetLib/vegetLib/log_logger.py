import logging

def log_make_logger(nameV):
   
    LOGGER = logging.getLogger(nameV)

    LOGGER.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
   
    myStreamTypicallySTDOUT = logging.StreamHandler()
    myStreamTypicallySTDOUT.setFormatter(formatter)
    LOGGER.addHandler(myStreamTypicallySTDOUT)

    fh = logging.FileHandler('./log/run.log')
    fh.setLevel(logging.INFO)
    fh.setFormatter(formatter)
    LOGGER.addHandler(fh)
    
    LOGGER.info("Logging Begins")
    return LOGGER

from inspect import currentframe

def log_get_line_number():
    cf = currentframe()
    return cf.f_back.f_lineno

