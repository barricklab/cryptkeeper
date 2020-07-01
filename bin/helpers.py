from time import time
import datetime

class Logger:
    '''Sends pretty colors to the console and also logs console to file'''
    def __init__(self):
        self.colors = {'normal': "\u001b[0m",
                  'warn': '\u001b[31m'}
        self.verbose = True
    @staticmethod
    def _send_to_log(text):
        pass
    def normal(self,text):
        if self.verbose:
            print(f"{self.colors['normal']}{text}{self.colors['normal']}")
        self._send_to_log(text)
    def warn(self, text):
        print(f"{self.colors['warn']}Warning: {text}{self.colors['normal']}")
        text = "Warning " + text
        self._send_to_log(text)



def timer_decorator(method):
    '''Function decorator to time functions'''
    def _timer(*args, **kw):
        start_time = datetime.datetime.now()
        result = method(*args, **kw)
        end_time = datetime.datetime.now()
        time_delta = end_time - start_time
        seconds = datetime.timedelta.total_seconds(time_delta)
        print(f'{method.__name__}  {seconds / 60} minutes')
        return result
    return _timer