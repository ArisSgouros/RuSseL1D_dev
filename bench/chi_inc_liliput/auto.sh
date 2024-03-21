RSL1D

merge_tool="../../tools/merge_prof.py"

python $merge_tool o.phi 0 1 # phi_mxa
python $merge_tool o.phi 0 2 # phi_mxb
python $merge_tool o.phi 0 7 # phi_tot
python $merge_tool o.field 0 1 # o.field_kd1
python $merge_tool o.field 0 2 # o.field_kd1_new
python $merge_tool o.field 0 3 # o.field_kd2
python $merge_tool o.field 0 4 # o.field_kd2_new
