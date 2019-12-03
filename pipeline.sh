#!/usr/bin/env bash
#
# Tiny and portable bash script as an executable command
#
# Usage:
#   pipeline.sh [--debug] [--cpus=<int>] [--bool] [--var=<str>] [<arg>...]
#   pipeline.sh --version
#   pipeline.sh -h|--help
#
# Options:
#   --debug           Debug mode
#   --cpus=<int>      Limit CPUs for multiprocessing
#   --bool            Option without an argument
#   --var=<str>       Option with an argument
#   --version         Print version
#   -h, --help        Print usage
#
# Arguments:
#   <arg>             Main arguments

set -ue

if [[ ${#} -ge 1 ]]; then
  for a in "${@}"; do
    [[ "${a}" = '--debug' ]] && set -x && break
  done
fi

COMMAND_PATH=$(realpath "${0}")
COMMAND_NAME=$(basename "${COMMAND_PATH}")
COMMAND_VERSION='v0.0.1'

case "${OSTYPE}" in
  darwin* )
    CPUS=$(sysctl -n hw.ncpu)
    ;;
  linux* )
    CPUS=$(grep -ce '^processor\s\+:' /proc/cpuinfo)
    ;;
  * )
    CPUS=1
    :
    ;;
esac
MAIN_ARGS=()
OPTION_BOOL=0
OPTION_VAR='NULL'

function print_version {
  echo "${COMMAND_NAME}: ${COMMAND_VERSION}"
}

function print_usage {
  sed -ne '1,2d; /^#/!q; s/^#$/# /; s/^# //p;' "${COMMAND_PATH}"
}

function abort {
  {
    if [[ ${#} -eq 0 ]]; then
      cat -
    else
      COMMAND_NAME=$(basename "${COMMAND_PATH}")
      echo "${COMMAND_NAME}: ${*}"
    fi
  } >&2
  exit 1
}

while [[ ${#} -ge 1 ]]; do
  case "${1}" in
    '--debug' )
      shift 1
      ;;
    '--bool' )
      OPTION_BOOL=1 && shift 1
      ;;
    '--var' )
      OPTION_VAR="${2}" && shift 2
      ;;
    --var=* )
      OPTION_VAR="${1#*\=}" && shift 1
      ;;
    '--cpus' )
      CPUS="${2}" && shift 2
      ;;
    --cpus=* )
      CPUS="${1#*\=}" && shift 1
      ;;
    '--version' )
      print_version && exit 0
      ;;
    '-h' | '--help' )
      print_usage && exit 0
      ;;
    -* )
      abort "invalid option: ${1}"
      ;;
    * )
      MAIN_ARGS+=("${1}") && shift 1
      ;;
  esac
done

[[ ${#MAIN_ARGS[@]} -gt 0 ]] \
  && ARGS_STR="${MAIN_ARGS[*]}" \
  || ARGS_STR='NULL'
printf "OSTYPE:      \t%s\n" "${OSTYPE}"
printf "BASH:        \t%s\n" "$(command -v bash)"
printf "COMMAND_PATH:\t%s\n" "${COMMAND_PATH}"
printf "PWD:         \t%s\n" "${PWD}"
printf "OPTION_BOOL: \t%d\n" "${OPTION_BOOL}"
printf "OPTION_VAR:  \t%s\n" "${OPTION_VAR}"
printf "CPUS:        \t%d\n" "${CPUS}"
printf "MAIN_ARGS:   \t%s\n" "${ARGS_STR}"
