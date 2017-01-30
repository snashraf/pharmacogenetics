#!/bin/sh
# Run a test for each testfile passed as an argument. If a testfile does not
# contain the expected test result, generate one.
set -e
trap 'rm -rf build/' INT
for TESTFILE; do
  printf 'Testfile %s\n' $TESTFILE
  for FORMAT in templates/*/; do
    printf '  Format %s\n' $FORMAT
    for TEMPLATE in ${FORMAT}*.tex.m4; do
      printf '    Template %s\n' $TEMPLATE
      m4 -DTEST_FILENAME=test.tex <$FORMAT/COMMANDS.m4 |
      (while read COMMAND; do
        printf '      Command %s\n' "$COMMAND"

        # Set up the testing directory.
        rm -rf build/
        mkdir build/
        cp support/* $TESTFILE build/
        cd build/
        sed -r '/^\s*<<<\s*$/{x;q}' \
          <${TESTFILE##*/} >test-setup.tex
        sed -rn '/^\s*<<<\s*$/,/^\s*>>>\s*$/{/^\s*(<<<|>>>)\s*$/!p}' \
          <${TESTFILE##*/} >test-input.md
        sed -n '/^\s*>>>\s*$/,${/^\s*>>>\s*$/!p}' \
          <${TESTFILE##*/} >test-expected.log
        m4 -DTEST_SETUP_FILENAME=test-setup.tex \
           -DTEST_INPUT_FILENAME=test-input.md <../$TEMPLATE >test.tex

        # Run the test, filter the output and concatenate adjacent lines.
        eval $COMMAND >/dev/null 2>&1
        sed -nr '/^\s*TEST INPUT BEGIN\s*$/,/^\s*TEST INPUT END\s*$/{
          /^\s*TEST INPUT (BEGIN|END)\s*$/!H
          /^\s*TEST INPUT END\s*$/{s/.*//;x;s/\n//g;p}
        }' <test.log >test-actual.log

        # If the testfile does not contain an expected outcome, use the current
        # outcome and update the testfile.
        if ! grep -q '^\s*>>>\s*$' <${TESTFILE##*/}; then
          cp test-actual.log test-expected.log
          (cat ${TESTFILE##*/} && printf '>>>\n'
            cat test-expected.log) >../$TESTFILE
          printf '      Added the expected test outcome to the testfile.\n'
        fi

        # Compare the expected outcome against the actual outcome.
        diff test-expected.log test-actual.log # ||
#          (sed -n '1,/^\s*>>>\s*$/p' <${TESTFILE##*/} && 
#            cat test-actual.log) >../$TESTFILE

        # Clean up the testing directory.
        cd ..
        rm -rf build/ # && break

      done) # && break
    done # && break
  done
done
