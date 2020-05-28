import logging

def log_make_logger(nameV):
   
    LOGGER = logging.getLogger(nameV)

    LOGGER.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
   
    myStreamTypicallySTDOUT = logging.StreamHandler()
    myStreamTypicallySTDOUT.setFormatter(formatter)
    LOGGER.addHandler(myStreamTypicallySTDOUT)
    
    LOGGER.info("Logging Begins")
    return LOGGER

from inspect import currentframe

def log_get_line_number():
    cf = currentframe()
    return cf.f_back.f_lineno

