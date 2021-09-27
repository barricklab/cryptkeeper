from time import time
import datetime
import logging


def timer_decorator(method):
    '''Function decorator to time functions'''
    def _timer(*args, **kw):
        start_time = datetime.datetime.now()
        result = method(*args, **kw)
        end_time = datetime.datetime.now()
        time_delta = end_time - start_time
        seconds = datetime.timedelta.total_seconds(time_delta)
        logging.DEBUG(f'{method.__name__}  {seconds / 60} minutes')
        return result
    return _timer