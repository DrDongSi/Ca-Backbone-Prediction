
--------------------------------------------------------------------------------
v1: is the first version that works.

*predict_set_cm.v1, predict_seq_cm.v1: from CMapPro, should be exactly same
with scratch server. (right now: use this version) 

predict_seq_cm.v2, predict_set_cm.v2: compile from cmappro (need to add lower_band,
upper_band to both model definition file and model file) [will be changed
soon]
now predict_seq_cm and predict_set_cm are compiled from cmappro. but unix
version is not updated yet. 

*predict_seq_sa: compile from acc-ensemble (model definition and model file don't need to be changed)

*predict_seq_ss: compile from sspro3_3d_comp (model definition file don't need
to be changed, model file need to be change(add three 0's, no use 3d, not use
acc, not use bin_num. )

*predict_band_cm.v1: predict contact map for a band. compiled form
crnn_known2
--------------------------------------------------------------------------------
