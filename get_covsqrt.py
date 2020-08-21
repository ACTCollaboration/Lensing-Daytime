import os,sys

# //// ACT mr3 for night time //// #
'''
for array in ['pa1','pa2','pa3']:
    for patch in ['deep56','boss']:
        cmd = "python actsims/bin/make_covsqrt.py v6.3.0_calibrated act_mr3 --season s15 --patch "+array+" --array "+array+" --covsqrt-kind multipow --overwrite --mask-version padded_v1 --nsims 3 --debug --delta-ell 100"
        os.system(cmd)
'''

# //// DR5 for daytime //// #
#qid = ['boss_d01','boss_d02']
#qid = ['boss_d03']
#qid = ['s16_d01','s16_d02']
qid = ['bndn_01','bndn_02']

for q in qid:
    # calculate cov_sqrt
    cmd = "python actsims/bin/make_covsqrt_dr5.py v6.3.0_calibrated dr5 --covsqrt-kind multipow --qid "+q+" --overwrite --fill-const --delta-ell 100 --nsims 0 --debug --calc_covsqrt"
    # run sim with pre-computed cov_sqrt
    #cmd = "python actsims/bin/make_covsqrt_dr5.py v6.3.0_calibrated dr5 --covsqrt-kind multipow --qid "+q+" --overwrite --fill-const --delta-ell 100 --nsims 3 --debug"
    os.system(cmd)


