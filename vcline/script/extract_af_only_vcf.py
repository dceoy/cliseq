#!/usr/bin/env python

import fileinput
import re
import sys


def main():
    search_regex = re.compile('\tPASS\t.*[\t;]AF=\\S')
    sub_regexes = [re.compile('\t\\S*;?(AF=[0-9e\\.\\+\\-]+)\\S*'), '\t\\1']
    for line in fileinput.input():
        if line.startswith('#'):
            sys.stdout.write(line)
        elif search_regex.search(line):
            sys.stdout.write(re.sub(*sub_regexes, line))


if __name__ == '__main__':
    main()
