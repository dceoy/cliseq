#!/usr/bin/env python

import fileinput
import os
import re
import sys


def main():
    search_regex = re.compile('@[A-Z]{1}')
    for line in fileinput.input():
        if not search_regex.match(line):
            v = line.split('\t')
            sys.stdout.write(
                '{0}\t{1}\t{2}'.format(v[0], (int(v[1]) - 1), int(v[2]))
                + os.linesep
            )
            sys.stdout.flush()


if __name__ == '__main__':
    main()
