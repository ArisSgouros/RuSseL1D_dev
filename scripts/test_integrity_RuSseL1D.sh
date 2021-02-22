#!/bin/bash
#
# The script works as follows:
# 1. loops over all tests declared in the TEST_LIST
# 2. copies the input files to the source
# 3. runs a simulation based on the input files
# 4. Compares the outputs specified in file LOG_LIST
#    located in each folder, with the outputs derived
#    from the simulation
#
#################################################################
EXEC_NAME="RuSseL1D"

# Beware! Set it to true only if you know what you know what you are doing!
replace_log=false
#replace_log=true

# Enable this and you will have the option to vimdiff each set of failed logs
inspect_log=true
#inspect_log=false

echo
if [ "$replace_log" = true ] ; then
   echo -e "Flag replace_log is set to true. You will have the option to \e[32mREPLACE\e[0m log files."
fi

#cd run

mv input.in.txt TEMP.input.in.txt

declare -a TEST_LIST=( \
   "test_integrity/0.ps_fg_sph_ham_it1/" \
   "test_integrity/1.ps_fg_sph_ham/" \
   "test_integrity/2.cafe222_ps_fg_sph_neumann_ham/" \
   "test_integrity/3.cafe222_ps_fg_planar_dirichlet_ham/" \
   "test_integrity/4.isolated/" \
   "test_integrity/5.isolated_planar/" \
   "test_integrity/6.ps_fg_sph_ham_nconst/" \
   "test_integrity/7.pe_f_SL_thomas2.0/" \
   "test_integrity/9.ps_fgg_HF_ham" \
   "test_integrity/10.test_with_tolis_wo_sqgrad" \
   "test_integrity/11.test_with_tolis_w_sqgrad" \
   "test_integrity/12.test_with_tolis_w_sqgrad_rect" \
   "test_integrity/13.ps_f_hybrid" \
)


echo 'Compiling..'
make cleantest >/dev/null
make >/dev/null

tests_failed=false
echo 'Initiating tests..'
echo '-----------------------------------------------------------------'
# Loop over all test cases
ii=0
for TEST_FOLDER in "${TEST_LIST[@]}"; do
   ii="$(( ii + 1 ))"
   echo
   echo -e "Testing folder: \e[36m$TEST_FOLDER\e[0m"

   # This reads the logs to be compared from the LOG_LISTS.txt
   # file in each folder.
   mapfile -t SOURCE_LIST < "$TEST_FOLDER"/SOURCE_LIST.txt

   for SOURCE in "${SOURCE_LIST[@]}"; do
      cp "$TEST_FOLDER"/"$SOURCE" "$SOURCE"
   done

   make cleanout

   echo "Running $EXEC_NAME.."

   if test ! -f "$EXEC_NAME"; then
      echo -e "Compilation \e[31mFAILED!\e[0m"
      echo -e "Current test is aborted.."
      continue;
   fi

   ./"$EXEC_NAME" > LOG.out.txt."$ii"
   #mpirun -np 2 ./"$EXEC_NAME" > LOG.out.txt."$ii"

   # This reads the logs to be compared from the LOG_LISTS.txt
   # file in each folder.
   mapfile -t LOG_LIST < "$TEST_FOLDER"/LOG_LIST.txt

   for LOG in "${LOG_LIST[@]}"; do

      if [ "$replace_log" = true ] ; then
         if test -f "$LOG"; then
            read -p "Replace log "$LOG"? (Y/n)" -n 1 -r
            echo
            if [[ $REPLY =~ ^[Yy]$ ]]; then
               cp "$LOG" "$TEST_FOLDER"/"$LOG"
               echo -e "LOG FILE \e[36m$LOG\e[0m has been \e[32mREPLACED!\e[0m"
            fi
            continue
         fi
      fi

      if test ! -f "$LOG"; then
         echo -e "LOG FILE: \e[36m$LOG\e[0m not found.. \e[93mTEST SKIPPED..\e[0m"
         tests_failed=true
         continue
      fi

      if test ! -f "$TEST_FOLDER"/"$LOG"; then
         echo -e "LOG FILE: \e[36m$TEST_FOLDER/$LOG\e[0m not found.. \e[30mTEST SKIPPED..\e[0m"
         tests_failed=true
         continue;
      fi

      OUTPUT="$(diff -q "$LOG" "$TEST_FOLDER"/"$LOG")"

      if [ -z "$OUTPUT" ]
      then
         echo -e "LOG FILE: \e[36m$LOG\e[0m .. \e[32mSUCCESS!\e[0m"
      else
         echo -e "LOG FILE: \e[36m$LOG\e[0m .. \e[31mFAILED!\e[0m"
         echo "${OUTPUT}"
         tests_failed=true

         if [ "$inspect_log" = true ] ; then
            read -p "Do you want to inspect the mismatched logs? (Y/n)" -n 1 -r
            echo
            if [[ $REPLY =~ ^[Yy]$ ]]; then
               vimdiff "$LOG" "$TEST_FOLDER"/"$LOG"
            fi
         fi
      fi
   done
done

if [ "$tests_failed" = true ] ; then
   echo -e '\e[36mNote: in case some tests are failed/skipped make sure the Makefile is set on DEBUG mode, MAKE_PRODUCTION_RUN=0\e[0m'
fi

mv TEMP.input.in.txt input.in.txt
