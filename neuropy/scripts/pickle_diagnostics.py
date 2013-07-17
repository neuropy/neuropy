"""Diagnose pickling errors of an obj, giving more detailed error messages.
See http://stackoverflow.com/a/7843281/2020363"""

import pickle
import StringIO
output = StringIO.StringIO()

class MyPickler (pickle.Pickler):
    """Not entirely sure how this works, but I think it recursively tries to pickle every
    attribute of obj, and prints out those that aren't picklable"""
    def save(self, obj):
        try:
            pickle.Pickler.save(self, obj)
        except Exception, e:
            # this print is optional, might want to disable at first if getting lots
            # of errors:
            print(obj)
            raise e

# diagnostic pickling:
import StringIO
output = StringIO.StringIO()
MyPickler(output).dump(obj)

# normal full speed pickling:
import cPickle
from io import BytesIO
cPickle.dump(obj, BytesIO(), cPickle.HIGHEST_PROTOCOL)
