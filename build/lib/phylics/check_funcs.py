# ==========================================================================
#                                  PhyliCS
# ==========================================================================
# This file is part of PhyliCS.
#
# TOOL is Free Software: you can redistribute it and/or modify it
# under the terms found in the LICENSE.rst file distributed
# together with this file.
#
# PhyliCS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
# ==========================================================================
# Author: Marilisa Montemurro <marilisa.montemurro@polito.it>
# ==========================================================================
# check_funcs.py: Check functions module
# ==========================================================================

import sys
import re
import argparse

#sys.tracebacklimit = 0

def check_valid_sample(nmin):
    class CheckValid(argparse.Action):
        def __call__(self, parser, args, values, option_string=None):
            if not nmin<=len(values):
                msg='argument "{f}" requires at least {nmin} arguments'.format(
                        f=self.dest,nmin=nmin)
                raise argparse.ArgumentTypeError(msg)
            for value in values:
                arguments = value.split(":")
                if not len(arguments) == 2:
                    msg='argument "{f}" must be formatted as follows: sample_name:input_dir'.format(f=self.dest)
                    raise argparse.ArgumentTypeError(msg)
            setattr(args, self.dest, values)
    return CheckValid

def check_valid_samples():
    class CheckValid(argparse.Action):
        def __call__(self, parser, args, values, option_string=None):
            for value in values:
                arguments = value.split(":")
                if not len(arguments) == 2:
                    msg='Wrong format: {}\nArgument {f} must be formatted as follows: sample_name:file_path'.format(value, f=self.dest)
                    raise argparse.ArgumentTypeError(msg)
            setattr(args, self.dest, values)
    return CheckValid

def check_valid_N_clust():
    class CheckValid(argparse.Action):
        def __call__(self, parser, args, value, option_string=None):
            if value <= 0:
                msg='argument "{f}" requires a number of clusters greater than 0'.format(f=self.dest)
                raise argparse.ArgumentTypeError(msg)
            setattr(args, self.dest, value)
    return CheckValid

def check_valid():
    class CheckValid(argparse.Action):
        def __call__(self, parser, args, values, option_string=None):
            if values[0] <= 0:
                msg='argument "{f}" requires a number of clusters greater than 0'.format(f=self.dest)
                raise argparse.ArgumentTypeError(msg)
            setattr(args, self.dest, values)
    return CheckValid

def check_valid_N(nmin):
    class CheckValid(argparse.Action):
        def __call__(self, parser, args, values, option_string=None):
            if not nmin<=len(values):
                msg='argument "{f}" requires at least {nmin} arguments'.format(f=self.dest, nmin=nmin)
                raise argparse.ArgumentTypeError('expected at least two arguments')
            setattr(args, self.dest, values)
    return CheckValid

def valid_interval():
    class CheckValid(argparse.Action):
        def __call__(self, parser, args, values, option_string=None):
            r = re.compile('^([0-9]+\.?[0-9]*)\-([0-9]+\.?[0-9]*)$')
            for v in values:
                msg = "Not a valid interval: '{0}'.".format(v)
                if r.match(v) == None:
                    raise argparse.ArgumentTypeError(msg)
                if len(v.split('-')) != 2:
                    raise argparse.ArgumentTypeError(msg)
                if float(v.split('-')[0]) >= float(v.split('-')[1]):
                    raise argparse.ArgumentTypeError(msg)
                setattr(args, self.dest, values)
    return CheckValid

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def is_interval(s):
    return len(s.split('-')) == 2

def check_valid_outliers():
    class CheckValid(argparse.Action):
        def __call__(self, parser, args, values, option_string=None):
            for v in values:
                fields = v.split('-')
                if len(fields) == 2:
                    if not is_number(fields[0]) or not is_number(fields[1]):
                        msg='[not a number] argument "{f}" requires a list of intervals or a single values. Intervals must be specified according to the following format: p1-p2, where p1 and p1 are both positive decimal numbers and p1 < p2. Single values must be specified as single positive decimal numbers'.format(f=self.dest)
                        raise argparse.ArgumentTypeError(msg)
                    if float(fields[0]) >= float(fields[1]):
                        msg='[p1>=p2] argument "{f}" requires a list of intervals or a single values. Intervals must be specified according to the following format: p1-p2, where p1 and p1 are both positive decimal numbers and p1 < p2. Single values must be specified as single positive decimal numbers'.format(f=self.dest)
                        raise argparse.ArgumentTypeError(msg)
                elif len(fields) == 1:
                    if not is_number(fields[0]):
                        msg='[not a number] argument "{f}" requires a list of intervals or a single values. Intervals must be specified according to the following format: p1-p2, where p1 and p1 are both positive decimal numbers and p1 < p2. Single values must be specified as single positive decimal numbers'.format(f=self.dest)
                        raise argparse.ArgumentTypeError(msg)
                else:
                    msg='[len={l}] argument "{f}" requires a list of intervals or a single values. Intervals must be specified according to the following format: p1-p2, where p1 and p1 are both positive decimal numbers and p1 < p2. Single values must be specified as single positive decimal numbers'.format(f=self.dest, l=len(fields))
                    raise argparse.ArgumentTypeError(msg)
            setattr(args, self.dest, values)
    return CheckValid

