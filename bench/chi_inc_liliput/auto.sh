RSL1D

merge_tool="../../tools/merge_prof.py"

#python $merge_tool o.phi 0 1 -shift -1 # phi_mxa
#python $merge_tool o.phi 0 2 # phi_mxb
python $merge_tool o.phi 0 7 -shift -1 -transpose 1 # phi_tot
python $merge_tool o.field 0 1 -transpose 1 # o.field_kd1
python $merge_tool o.field 0 2 -transpose 1 # o.field_kd1_new
python $merge_tool o.field 0 3 -transpose 1 # o.field_kd2
python $merge_tool o.field 0 4 -transpose 1 # o.field_kd2_new

paste phi_tot.csv wa_kd1.csv wa_kd2.csv wa_new_kd1.csv wa_new_kd2.csv > all.paste.csv

cp *.csv /media/share/chi_inc_liliput/.
