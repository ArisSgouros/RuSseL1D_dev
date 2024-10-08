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

root_path=$(pwd)"/"

run_path=$root_path"/run/"
test_path=$root_path"/test_integrity/"
exec_path=$run_path"RSL1D"

#exit
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

mv $run_path"/in.input" $run_path"/TEMP.in.input"

declare -a TEST_LIST=( \
   $test_path"/16.flory_huggins_graft" \
   $test_path"/15.flory_huggins_short" \
   $test_path"/14.flory_huggins" \
   $test_path"/13.8_w_matrixB/" \
   $test_path"/0.ps_fg_sph_ham_it1/" \
   $test_path"/1.ps_fg_sph_ham/" \
   $test_path"/2.cafe222_ps_fg_sph_neumann_ham/" \
   $test_path"/3.cafe222_ps_fg_planar_dirichlet_ham/" \
   $test_path"/4.isolated/" \
   $test_path"/5.isolated_planar/" \
   $test_path"/6.ps_fg_sph_ham_nconst/" \
   $test_path"/7.pe_f_SL_thomas2.0/" \
   $test_path"/8.ps_fgg_HF_ham" \
   $test_path"/9.test_with_tolis_wo_sqgrad" \
   $test_path"/10.test_with_tolis_w_sqgrad" \
   $test_path"/11.test_with_tolis_w_sqgrad_rect" \
   $test_path"/12.ps_f_hybrid" \
)

: '


'

echo 'Compiling..'
make cleantest >/dev/null
make >/dev/null

cd $run_path

tests_failed=false
echo 'Initiating tests..'
echo '-----------------------------------------------------------------'
# Loop over all test cases
ii=0
for TEST_FOLDER in "${TEST_LIST[@]}"; do
   ii="$(( ii + 1 ))"
   echo
   echo -e "Testing folder: \e[36m$TEST_FOLDER\e[0m"

   # This reads the logs to be compared from the LOG_LISTS
   # file in each folder.
   mapfile -t SOURCE_LIST < "$TEST_FOLDER"/SOURCE_LIST

   for SOURCE in "${SOURCE_LIST[@]}"; do
      cp "$TEST_FOLDER"/"$SOURCE" "$SOURCE"
   done

   cd $root_path
   make cleanout
   cd $run_path

   echo "Running $exec_path.."

   if test ! -f "$exec_path"; then
      echo -e "Compilation \e[31mFAILED!\e[0m"
      echo -e "Current test is aborted.."
      continue;
   fi

   "$exec_path" > LOG."$ii"
   #mpirun -np 2 ./"$exec_path" > LOG."$ii"

   # This reads the logs to be compared from the LOG_LISTS
   # file in each folder.
   mapfile -t LOG_LIST < "$TEST_FOLDER"/LOG_LIST

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

cd $run_path
mv TEMP.in.input in.input
