#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (C) 2016 INRA
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

__author__ = 'Maria Bernard - Sigenae Team'
__copyright__ = 'Copyright (C) 2016 INRA'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'maria.bernard@inra.fr'

import sys
import argparse

from stacks.stacks_summary import *

### MAIN test-data/test_results/stacks_outputs/summary.html
def __main__():
    parser = argparse.ArgumentParser()
    parser.add_argument('--stacks-prog', choices=["ustacks","pstacks","cstacks","sstacks","genotypes","populations","ref_map.pl","denovo_map.pl"], required=True)
    parser.add_argument('--res-dir', required=True)
    parser.add_argument('--logfile', default=None, required=False)
    parser.add_argument('--pop-map', default=None)
    parser.add_argument('--summary', required=True)
    options = parser.parse_args()
    print "Command : " + " ".join(sys.argv)
    summarize_results( options.stacks_prog, options.res_dir, options.logfile, options.pop_map, options.summary )


if __name__ == '__main__':
    __main__()
