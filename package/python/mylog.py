import sys
import logging

class StreamToLogger(object):
    """
    Fake file-like stream object that redirects writes to a logger instance.
    Ref:
    https://stackoverflow.com/questions/19425736/how-to-redirect-stdout-and-stderr-to-logger-in-python
    """
    def __init__(self, logger, level):
       self.logger = logger
       self.level = level
       self.linebuf = ''

    def write(self, buf):
       for line in buf.rstrip().splitlines():
          self.logger.log(self.level, line.rstrip())

    def flush(self):
        pass

def set_file_logger(fnm:str,
                    fmode:str = 'a',
                    name:str = 'sa2_pp',
                    log_level: int = logging.DEBUG) -> logging.Logger:
    logger = logging.getLogger(name)
    logger.setLevel(log_level)
    fh = logging.FileHandler(filename = fnm,
                             mode = fmode)
    fm = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fh.setFormatter(fm)
    logger.addHandler(fh)
    return logger

def handle_exception(logger, exc_type, exc_value, exc_traceback):
    import traceback
    if issubclass(exc_type, KeyboardInterrupt):
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return
    logger.error(''.join(["Uncaught exception: ",
                         *traceback.format_exception(
                             exc_type, exc_value, exc_traceback)
                         ]))
