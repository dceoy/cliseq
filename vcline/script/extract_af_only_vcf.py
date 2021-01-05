#!/usr/bin/env python

import fileinput
import re
import sys


def main():
    search_regex = re.compile('\tPASS\t.*[\t;]AF=\\S')
    sub_regex = re.compile('\t\\S*;?(AF=[0-9e\\.\\+\\-]+)\\S*')
    for line in fileinput.input():
        if line.startswith('#'):
            sys.stdout.write(line)
        elif search_regex.search(line):
            sys.stdout.write(re.sub(sub_regex, '\t\\1', line))


if __name__ == '__main__':
    main()
