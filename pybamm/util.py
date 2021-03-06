#
# Utility classes for PyBaMM
#
# The code in this file is adapted from Pints
# (see https://github.com/pints-team/pints)
#
from __future__ import absolute_import, division
from __future__ import print_function, unicode_literals

import cProfile
import importlib
import os
import pstats
import sys
import timeit


class Timer(object):
    """
    Provides accurate timing.

    Example
    -------
    timer = pybamm.Timer()
    print(timer.format(timer.time()))

    """

    def __init__(self):
        self._start = timeit.default_timer()

    def format(self, time=None):
        """
        Formats a (non-integer) number of seconds, returns a string like
        "5 weeks, 3 days, 1 hour, 4 minutes, 9 seconds", or "0.0019 seconds".

        Arguments
        ---------
        time : float, optional
            The time to be formatted.

        Returns
        -------
        string
            The string representation of ``time`` in human-readable form.
        """
        if time is None:
            time = self.time()
        if time < 1e-2:
            return str(time) + " seconds"
        elif time < 60:
            return str(round(time, 2)) + " seconds"
        output = []
        time = int(round(time))
        units = [(604800, "week"), (86400, "day"), (3600, "hour"), (60, "minute")]
        for k, name in units:
            f = time // k
            if f > 0 or output:
                output.append(str(f) + " " + (name if f == 1 else name + "s"))
            time -= f * k
        output.append("1 second" if time == 1 else str(time) + " seconds")
        return ", ".join(output)

    def reset(self):
        """
        Resets this timer's start time.
        """
        self._start = timeit.default_timer()

    def time(self):
        """
        Returns the time (float, in seconds) since this timer was created,
        or since meth:`reset()` was last called.
        """
        return timeit.default_timer() - self._start


def profile(code, sort="cumulative", num=30):
    """Common-use for cProfile"""
    cProfile.run(code)
    stats = pstats.Stats()
    stats.sort_stats(sort)
    return stats


def load_function(filename):
    """
    Load a python function from a file "function_name.py" called "function_name".
    The filename might either be an absolute path, in which case that specific file will
    be used, or the file will be searched for relative to PyBaMM root.

    Arguments
    ---------
    filename : str
        The name of the file containing the function of the same name.

    Returns
    -------
    function
        The python function loaded from the file.
    """

    if not filename.endswith('.py'):
        raise ValueError("Expected filename.py, but got {}".format(filename))

    # If it's an absolute path, find that exact file
    if os.path.isabs(filename):
        if not os.path.isfile(filename):
            raise ValueError(
                "{} is an absolute path, but the file is not found".format(filename))

        valid_filename = filename

    # Else, search in the whole PyBaMM directory for matches
    else:
        search_path = os.path.commonpath([pth for pth in sys.path if 'PyBaMM' in pth])
        head, tail = os.path.split(filename)

        matching_files = []

        for root, _, files in os.walk(search_path):
            for file in files:
                if file == tail:
                    full_path = os.path.join(root, file)
                    if full_path.endswith(filename):
                        matching_files.append(full_path)

        if len(matching_files) == 0:
            raise ValueError(
                "{} cannot be found in the PyBaMM directory".format(filename))
        elif len(matching_files) > 1:
            raise ValueError(
                "{} found multiple times in the PyBaMM directory".format(filename))

        valid_filename = matching_files[0]

    # Now: we have some /path/to/valid/filename.py
    # Add "/path/to/vaid" to the python path, and load the module "filename".
    # Then, check "filename" module contains "filename" function.  If it does, return
    # that function object, or raise an exception

    valid_path, valid_leaf = os.path.split(valid_filename)
    sys.path.append(valid_path)

    # Load the module, which must be the leaf of filename, minus the .py extension
    valid_module = valid_leaf.replace('.py', '')
    module_object = importlib.import_module(valid_module)

    # Check that a function of the same name exists in the loaded module
    if valid_module not in dir(module_object):
        raise ValueError(
            "No function {} found in module {}".format(valid_module, valid_module))

    return getattr(module_object, valid_module)
