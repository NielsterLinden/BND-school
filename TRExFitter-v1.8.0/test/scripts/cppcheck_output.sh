#!/bin/bash

EXISTING_WARNINGS=0
EXISTING_STYLE=0
EXISTING_PERFORMANCE=0
EXISTING_PORTABILITY=0
EXISTING_INFORMATION=1

if [[ $(cat cppcheck.txt | wc -l) -gt 0 ]]; then
  # if the file contains any lines
  ERRORS=$(cat cppcheck.txt | grep ": error:" | wc -l)
  WARNINGS=$(cat cppcheck.txt | grep ": warning:" | wc -l)
  STYLE=$(cat cppcheck.txt | grep ": style:" | wc -l)
  PERFORMANCE=$(cat cppcheck.txt | grep ": performance:" | wc -l)
  PORTABILITY=$(cat cppcheck.txt | grep ": portability:" | wc -l)
  INFORMATION=$(cat cppcheck.txt | grep ": information:" | wc -l)
  echo "----------------------------------------------------"
  echo "cppcheck found the following amount of notifications:"
  echo "error:       " $ERRORS
  echo "warning:     " $WARNINGS
  echo "style:       " $STYLE
  echo "performance: " $PERFORMANCE
  echo "portability: " $PORTABILITY
  echo "information: " $INFORMATION

  # fail if any errors are found
  if [[ $ERRORS -gt 0 ]]; then
    echo $ERRORS "error(s) found"
    exit 1
  fi;

  # fail if number of warnings is increased, notify to update script if number is decreased
  if [[ $WARNINGS -gt $EXISTING_WARNINGS ]]; then
    echo "number of warnings increased, please fix"
    exit 1
  elif [[ $WARNINGS -lt $EXISTING_WARNINGS ]]; then
    echo "number of warnings decreased, please update test/scripts/cppcheck_output.sh"
    exit 1
  fi;

  # fail if number of style comments is increased, notify to update script if number is decreased
  if [[ $STYLE -gt $EXISTING_STYLE ]]; then
    echo "number of style comments increased, please fix"
    exit 1
  elif [[ $STYLE -lt $EXISTING_STYLE ]]; then
    echo "number of style comments decreased, please update test/scripts/cppcheck_output.sh"
    exit 1
  fi;

  # fail if number of performance comments is increased, notify to update script if number is decreased
  if [[ $PERFORMANCE -gt $EXISTING_PERFORMANCE ]]; then
    echo "number of performance comments increased, please fix"
    exit 1
  elif [[ $PERFORMANCE -lt $EXISTING_PERFORMANCE ]]; then
    echo "number of performance comments decreased, please update test/scripts/cppcheck_output.sh"
    exit 1
  fi;

  # fail if number of portability comments is increased, notify to update script if number is decreased
  if [[ $PORTABILITY -gt $EXISTING_PORTABILITY ]]; then
    echo "number of portability comments increased, please fix"
    exit 1
  elif [[ $PORTABILITY -lt $EXISTING_PORTABILITY ]]; then
    echo "number of portability comments decreased, please update test/scripts/cppcheck_output.sh"
    exit 1
  fi;

  # fail if number of  comments is increased, notify to update script if number is decreased
  if [[ $INFORMATION -gt $EXISTING_INFORMATION ]]; then
    echo "number of information comments increased, please fix"
    exit 1
  elif [[ $INFORMATION -lt $EXISTING_INFORMATION ]]; then
    echo "number of information comments decreased, please update test/scripts/cppcheck_output.sh"
    exit 1
  fi;

  # all checks passed successfully
  exit 0
fi;
