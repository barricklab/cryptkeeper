from time import time
import datetime
import shelve
import os
import inspect
import hashlib
from functools import wraps
from pathlib import Path

def timer(func):
    """ A timer decorator """
    def wrapper(*args, **kwargs):
        start = time()
        result = func(*args, **kwargs)
        end = time()
        print(f"{func.__name__} took {end-start} seconds.")
        return result
    return wrapper


def persistant_cache(filepath=None):
    """ A persistant cache decorator """
    def decorator(function):
        @wraps(function)
        def wrapper(*args, **kwargs):
            nonlocal filepath
            if filepath is None:
                filepath = "cache.db"
            with shelve.open(filepath) as cache:
                key = hashlib.md5(str("func + args + kwargs").encode('utf-8')).hexdigest()
                if key not in cache:
                        cache[key] = function(*args, **kwargs)
                return cache[key]
        return wrapper
    return decorator
