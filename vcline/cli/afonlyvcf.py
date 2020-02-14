#!/usr/bin/env python

import fileinput
import re
import sys


def main():
    search_regex = re.compile('\tPASS\t.*[\t;]AF=[^;]')
    sub_regexes = [
        re.compile('(\t[^\t]*;|\t)(AF=[0-9]*\\.[e0-9+-]*)[^\t\r\n]*'), '\t\\2'
    ]
    for line in fileinput.input():
        if line.startswith('#'):
            sys.stdout.write(line)
            sys.stdout.flush()
        elif search_regex.search(line):
            sys.stdout.write(re.sub(*sub_regexes, line))
            sys.stdout.flush()


if __name__ == '__main__':
    main()
