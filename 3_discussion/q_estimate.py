import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import binned_statistic
import batman
import emcee
import corner
import os, sys, time
from scipy.optimize import minimize


def lnlike(theta,                       orbit_number_1, t1, uncertainty_1,
                                        orbit_number_2, t2, uncertainty_2,
                                        orbit_number_3, t3, uncertainty_3,
                                        orbit_number_4, t4, uncertainty_4,
                                        orbit_number_5, t5, uncertainty_5,
                                        orbit_number_6, t6, uncertainty_6,
                                        orbit_number_7, t7, uncertainty_7,
                                        orbit_number_8, t8, uncertainty_8,
                                        orbit_number_9, t9, uncertainty_9,
                                        orbit_number_10, t10, uncertainty_10,
                                        orbit_number_11, t11, uncertainty_11,
                                        orbit_number_12, t12, uncertainty_12,
                                        orbit_number_13, t13, uncertainty_13,
                                        orbit_number_14, t14, uncertainty_14,
                                        orbit_number_15, t15, uncertainty_15,
                                        orbit_number_16, t16, uncertainty_16,
                                        orbit_number_17, t17, uncertainty_17,
                                        orbit_number_18, t18, uncertainty_18,
                                        orbit_number_19, t19, uncertainty_19,
                                        orbit_number_20, t20, uncertainty_20,
                                        orbit_number_21, t21, uncertainty_21,
                                        orbit_number_22, t22, uncertainty_22,
                                        orbit_number_23, t23, uncertainty_23,
                                        orbit_number_24, t24, uncertainty_24,
                                        orbit_number_25, t25, uncertainty_25,
                                        orbit_number_26, t26, uncertainty_26,
                                        orbit_number_27, t27, uncertainty_27,
                                        orbit_number_28, t28, uncertainty_28,
                                        orbit_number_29, t29, uncertainty_29,
                                        orbit_number_30, t30, uncertainty_30,
                                        factor_1,
                                        factor_2,
                                        factor_3,
                                        factor_4,
                                        factor_5,
                                        factor_6,
                                        factor_7,
                                        factor_8,
                                        factor_9,
                                        factor_10,
                                        factor_11,
                                        factor_12,
                                        factor_13,
                                        factor_14,
                                        factor_15,
                                        factor_16,
                                        factor_17,
                                        factor_18,
                                        factor_19,
                                        factor_20,
                                        factor_21,
                                        factor_22,
                                        factor_23,
                                        factor_24,
                                        factor_25,
                                        factor_26,
                                        factor_27,
                                        factor_28,
                                        factor_29,
                                        factor_30
                                        ):
  period_1, t0_1, period_2, t0_2, period_3, t0_3, \
  period_4, t0_4, period_5, t0_5, period_6, t0_6, \
  period_7, t0_7, period_8, t0_8, period_9, t0_9, \
  period_10, t0_10, period_11, t0_11, period_12, t0_12, \
  period_13, t0_13, period_14, t0_14, period_15, t0_15, \
  period_16, t0_16, period_17, t0_17, period_18, t0_18, \
  period_19, t0_19, period_20, t0_20, period_21, t0_21, \
  period_22, t0_22, period_23, t0_23, period_24, t0_24, \
  period_25, t0_25, period_26, t0_26, period_27, t0_27, \
  period_28, t0_28, period_29, t0_29, period_30, t0_30, q  = theta   


  model_1 = t0_1+ orbit_number_1*period_1 + 1/q * factor_1 * orbit_number_1**2 /2
  model_2 = t0_2+ orbit_number_2*period_2 + 1/q * factor_2 * orbit_number_2**2 /2
  model_3 = t0_3+ orbit_number_3*period_3 + 1/q * factor_3 * orbit_number_3**2 /2
  model_4 = t0_4+ orbit_number_4*period_4 + 1/q * factor_4 * orbit_number_4**2 /2
  model_5 = t0_5+ orbit_number_5*period_5 + 1/q * factor_5 * orbit_number_5**2 /2
  model_6 = t0_6+ orbit_number_6*period_6 + 1/q * factor_6 * orbit_number_6**2 /2
  model_7 = t0_7+ orbit_number_7*period_7 + 1/q * factor_7* orbit_number_7**2 /2
  model_8 = t0_8+ orbit_number_8*period_8 + 1/q * factor_8* orbit_number_8**2 /2
  model_9 = t0_9+ orbit_number_9*period_9 + 1/q * factor_9 * orbit_number_9**2 /2
  model_10 = t0_10+ orbit_number_10*period_10 + 1/q * factor_10 * orbit_number_10**2 /2
  model_11 = t0_11+ orbit_number_11*period_11 + 1/q * factor_11 * orbit_number_11**2 /2
  model_12 = t0_12+ orbit_number_12*period_12 + 1/q * factor_12* orbit_number_12**2 /2
  model_13 = t0_13+ orbit_number_13*period_13 + 1/q * factor_13* orbit_number_13**2 /2
  model_14 = t0_14+ orbit_number_14*period_14 + 1/q * factor_14* orbit_number_14**2 /2
  model_15 = t0_15+ orbit_number_15*period_15 + 1/q * factor_15* orbit_number_15**2 /2
  model_16 = t0_16+ orbit_number_16*period_16 + 1/q * factor_16* orbit_number_16**2 /2
  model_17 = t0_17+ orbit_number_17*period_17 + 1/q * factor_17* orbit_number_17**2 /2
  model_18 = t0_18+ orbit_number_18*period_18 + 1/q * factor_18* orbit_number_18**2 /2
  model_19 = t0_19+ orbit_number_19*period_19 + 1/q * factor_19* orbit_number_19**2 /2
  model_20 = t0_20+ orbit_number_20*period_20 + 1/q * factor_20* orbit_number_20**2 /2
  model_21 = t0_21+ orbit_number_21*period_21 + 1/q * factor_21* orbit_number_21**2 /2
  model_22 = t0_22+ orbit_number_23*period_22 + 1/q * factor_22* orbit_number_22**2 /2
  model_23 = t0_23+ orbit_number_23*period_23 + 1/q * factor_23* orbit_number_23**2 /2
  model_24 = t0_24+ orbit_number_24*period_24 + 1/q * factor_24* orbit_number_24**2 /2
  model_25 = t0_25+ orbit_number_25*period_25 + 1/q * factor_25* orbit_number_25**2 /2
  model_26 = t0_26+ orbit_number_26*period_26 + 1/q * factor_26* orbit_number_26**2 /2
  model_27 = t0_27+ orbit_number_27*period_27 + 1/q * factor_27* orbit_number_27**2 /2
  model_28 = t0_28+ orbit_number_28*period_28 + 1/q * factor_28* orbit_number_28**2 /2
  model_29 = t0_29+ orbit_number_29*period_29 + 1/q * factor_29* orbit_number_29**2 /2
  model_30 = t0_30+ orbit_number_30*period_30 + 1/q * factor_30* orbit_number_30**2 /2

  diff_1 = (t1 - model_1)**2/uncertainty_1**2
  diff_2 = (t2 - model_2)**2/uncertainty_2**2
  diff_3 = (t3 - model_3)**2/uncertainty_3**2
  diff_4 = (t4 - model_4)**2/uncertainty_4**2
  diff_5 = (t5 - model_5)**2/uncertainty_5**2
  diff_6 = (t6 - model_6)**2/uncertainty_6**2
  diff_7 = (t7 - model_7)**2/uncertainty_7**2
  diff_8 = (t8 - model_8)**2/uncertainty_8**2
  diff_9 = (t9 - model_9)**2/uncertainty_9**2
  diff_10 = (t10 - model_10)**2/uncertainty_10**2
  diff_11 = (t11 - model_11)**2/uncertainty_11**2
  diff_12 = (t12 - model_12)**2/uncertainty_12**2
  diff_13 = (t13 - model_13)**2/uncertainty_13**2
  diff_14 = (t14 - model_14)**2/uncertainty_14**2
  diff_15 = (t15 - model_15)**2/uncertainty_15**2
  diff_16 = (t16 - model_16)**2/uncertainty_16**2
  diff_17 = (t17 - model_17)**2/uncertainty_17**2
  diff_18 = (t18 - model_18)**2/uncertainty_18**2
  diff_19 = (t19 - model_19)**2/uncertainty_19**2
  diff_20 = (t20 - model_20)**2/uncertainty_20**2
  diff_21 = (t21 - model_21)**2/uncertainty_21**2
  diff_22 = (t22 - model_22)**2/uncertainty_22**2
  diff_23 = (t23 - model_23)**2/uncertainty_23**2
  diff_24 = (t24 - model_24)**2/uncertainty_24**2
  diff_25 = (t25 - model_25)**2/uncertainty_25**2
  diff_26 = (t26 - model_26)**2/uncertainty_26**2
  diff_27 = (t27 - model_27)**2/uncertainty_27**2
  diff_28 = (t28 - model_28)**2/uncertainty_28**2
  diff_29 = (t29 - model_29)**2/uncertainty_29**2
  diff_30 = (t30 - model_30)**2/uncertainty_30**2

  return -(np.sum(diff_1) + np.sum(diff_2)  + np.sum(diff_3)  + np.sum(diff_4)  + np.sum(diff_5) \
     + np.sum(diff_6) + np.sum(diff_7) + np.sum(diff_8) + np.sum(diff_9) + np.sum(diff_10) + np.sum(diff_11) \
      + np.sum(diff_12) + np.sum(diff_13) + np.sum(diff_14) + np.sum(diff_15) + np.sum(diff_16) + np.sum(diff_17) \
       + np.sum(diff_18) + np.sum(diff_19) + np.sum(diff_20) + np.sum(diff_21) + np.sum(diff_22) + np.sum(diff_23) \
        + np.sum(diff_24) + np.sum(diff_25) + np.sum(diff_26) + np.sum(diff_27) + np.sum(diff_28) + np.sum(diff_29) \
         + np.sum(diff_30))



def fold_time(t, t0, period):
    
    # returns folded time in between -period/2 and +period/2
    # as well as an orbit number, assigning t0 = orbit zero

    orbit_number, t_fold = np.divmod(t-(t0-0.5*period), period)
    t_fold = t_fold - 0.5*period
    
    return orbit_number.astype(np.int)


def main():
    planet_1 = 'HATS-70'
    planet_2 = 'WASP-019'
    planet_3 = 'WASP-018'
    planet_4 = 'WASP-103'
    planet_5 = 'HATS-18'
    planet_6 = 'WASP-077'
    planet_7 = 'WASP-043'
    planet_8 = 'WASP-012'
    planet_9 = 'KELT-16'
    planet_10 = 'TrES-5'
    planet_11 = 'OGLE-TR-056'
    planet_12 = 'WASP-121'
    planet_13 = 'WASP-033'
    planet_14 = 'HATS-24'
    planet_15 = 'Qatar-2'
    planet_16 = 'TrES-3'
    planet_17 = 'CoRoT-02'
    planet_18 = 'HATS-52'
    planet_19 = 'WASP-004'
    planet_20 = 'WASP-173'
    planet_21 = 'CoRoT-01'
    planet_22 = 'HAT-P-36'
    planet_23 = 'WASP-135'
    planet_24 = 'HAT-P-23'
    planet_25 = 'KELT-09'
    planet_26 = 'WASP-046'
    planet_27 = 'Qatar-4'
    planet_28 = 'WASP-122'
    planet_29 = 'WASP-036'
    planet_30 = 'HATS-23'

    path_1 = f'/Users/kate/Documents/research/paper/1_database/reviewed_planets_clean/{planet_1}-mid-times/'
    path_2 = f'/Users/kate/Documents/research/paper/1_database/reviewed_planets_clean/{planet_2}-mid-times/'
    path_3 = f'/Users/kate/Documents/research/paper/1_database/reviewed_planets_clean/{planet_3}-mid-times/'
    path_4 = f'/Users/kate/Documents/research/paper/1_database/reviewed_planets_clean/{planet_4}-mid-times/'
    path_5 = f'/Users/kate/Documents/research/paper/1_database/reviewed_planets_clean/{planet_5}-mid-times/'
    path_6 = f'/Users/kate/Documents/research/paper/1_database/reviewed_planets_clean/{planet_6}-mid-times/'
    path_7 = f'/Users/kate/Documents/research/paper/1_database/reviewed_planets_clean/{planet_7}-mid-times/'
    path_8 = f'/Users/kate/Documents/research/paper/1_database/reviewed_planets_clean/{planet_8}-mid-times/'
    path_9 = f'/Users/kate/Documents/research/paper/1_database/reviewed_planets_clean/{planet_9}-mid-times/'
    path_10 = f'/Users/kate/Documents/research/paper/1_database/reviewed_planets_clean/{planet_10}-mid-times/'
    path_11 = f'/Users/kate/Documents/research/paper/1_database/reviewed_planets_clean/{planet_11}-mid-times/'
    path_12 = f'/Users/kate/Documents/research/paper/1_database/reviewed_planets_clean/{planet_12}-mid-times/'
    path_13 = f'/Users/kate/Documents/research/paper/1_database/reviewed_planets_clean/{planet_13}-mid-times/'
    path_14 = f'/Users/kate/Documents/research/paper/1_database/reviewed_planets_clean/{planet_14}-mid-times/'
    path_15 = f'/Users/kate/Documents/research/paper/1_database/reviewed_planets_clean/{planet_15}-mid-times/'
    path_16 = f'/Users/kate/Documents/research/paper/1_database/reviewed_planets_clean/{planet_16}-mid-times/'
    path_17 = f'/Users/kate/Documents/research/paper/1_database/reviewed_planets_clean/{planet_17}-mid-times/'
    path_18 = f'/Users/kate/Documents/research/paper/1_database/reviewed_planets_clean/{planet_18}-mid-times/'
    path_19 = f'/Users/kate/Documents/research/paper/1_database/reviewed_planets_clean/{planet_19}-mid-times/'
    path_20 = f'/Users/kate/Documents/research/paper/1_database/reviewed_planets_clean/{planet_20}-mid-times/'
    path_21 = f'/Users/kate/Documents/research/paper/1_database/reviewed_planets_clean/{planet_21}-mid-times/'
    path_22 = f'/Users/kate/Documents/research/paper/1_database/reviewed_planets_clean/{planet_22}-mid-times/'
    path_23 = f'/Users/kate/Documents/research/paper/1_database/reviewed_planets_clean/{planet_23}-mid-times/'
    path_24 = f'/Users/kate/Documents/research/paper/1_database/reviewed_planets_clean/{planet_24}-mid-times/'
    path_25 = f'/Users/kate/Documents/research/paper/1_database/reviewed_planets_clean/{planet_25}-mid-times/'
    path_26 = f'/Users/kate/Documents/research/paper/1_database/reviewed_planets_clean/{planet_26}-mid-times/'
    path_27 = f'/Users/kate/Documents/research/paper/1_database/reviewed_planets_clean/{planet_27}-mid-times/'
    path_28 = f'/Users/kate/Documents/research/paper/1_database/reviewed_planets_clean/{planet_28}-mid-times/'
    path_29 = f'/Users/kate/Documents/research/paper/1_database/reviewed_planets_clean/{planet_29}-mid-times/'
    path_30 = f'/Users/kate/Documents/research/paper/1_database/reviewed_planets_clean/{planet_30}-mid-times/'

    path2ratio = '/Users/kate/Documents/research/paper/3_tables/1_target_list_tic_tmag_snr_mr_ratio.csv'
    df = pd.read_csv(path2ratio)

    factor_1 = df[df['System'] == planet_1]['(m_p/m_s)*(R_p/a)^5'].iloc[0] * (-27 * np.pi /2)
    factor_2 = df[df['System'] == planet_2]['(m_p/m_s)*(R_p/a)^5'].iloc[0] * (-27 * np.pi /2)
    factor_3 = df[df['System'] == planet_3]['(m_p/m_s)*(R_p/a)^5'].iloc[0] * (-27 * np.pi /2)
    factor_4 = df[df['System'] == planet_4]['(m_p/m_s)*(R_p/a)^5'].iloc[0] * (-27 * np.pi /2)
    factor_5 = df[df['System'] == planet_5]['(m_p/m_s)*(R_p/a)^5'].iloc[0] * (-27 * np.pi /2)
    factor_6 = df[df['System'] == planet_6]['(m_p/m_s)*(R_p/a)^5'].iloc[0] * (-27 * np.pi /2)
    factor_7 = df[df['System'] == planet_7]['(m_p/m_s)*(R_p/a)^5'].iloc[0] * (-27 * np.pi /2)
    factor_8 = df[df['System'] == planet_8]['(m_p/m_s)*(R_p/a)^5'].iloc[0] * (-27 * np.pi /2)
    factor_9 = df[df['System'] == planet_9]['(m_p/m_s)*(R_p/a)^5'].iloc[0] * (-27 * np.pi /2)
    factor_10 = df[df['System'] == planet_10]['(m_p/m_s)*(R_p/a)^5'].iloc[0] * (-27 * np.pi /2)
    factor_11 = df[df['System'] == planet_11]['(m_p/m_s)*(R_p/a)^5'].iloc[0] * (-27 * np.pi /2)
    factor_12 = df[df['System'] == planet_12]['(m_p/m_s)*(R_p/a)^5'].iloc[0] * (-27 * np.pi /2)
    factor_13 = df[df['System'] == planet_13]['(m_p/m_s)*(R_p/a)^5'].iloc[0] * (-27 * np.pi /2)
    factor_14 = df[df['System'] == planet_14]['(m_p/m_s)*(R_p/a)^5'].iloc[0] * (-27 * np.pi /2)
    factor_15 = df[df['System'] == planet_15]['(m_p/m_s)*(R_p/a)^5'].iloc[0] * (-27 * np.pi /2)
    factor_16 = df[df['System'] == planet_16]['(m_p/m_s)*(R_p/a)^5'].iloc[0] * (-27 * np.pi /2)
    factor_17 = df[df['System'] == planet_17]['(m_p/m_s)*(R_p/a)^5'].iloc[0] * (-27 * np.pi /2)
    factor_18 = df[df['System'] == planet_18]['(m_p/m_s)*(R_p/a)^5'].iloc[0] * (-27 * np.pi /2)
    factor_19 = df[df['System'] == planet_19]['(m_p/m_s)*(R_p/a)^5'].iloc[0] * (-27 * np.pi /2)
    factor_20 = df[df['System'] == planet_20]['(m_p/m_s)*(R_p/a)^5'].iloc[0] * (-27 * np.pi /2)
    factor_21 = df[df['System'] == planet_21]['(m_p/m_s)*(R_p/a)^5'].iloc[0] * (-27 * np.pi /2)
    factor_22 = df[df['System'] == planet_22]['(m_p/m_s)*(R_p/a)^5'].iloc[0] * (-27 * np.pi /2)
    factor_23 = df[df['System'] == planet_23]['(m_p/m_s)*(R_p/a)^5'].iloc[0] * (-27 * np.pi /2)
    factor_24 = df[df['System'] == planet_24]['(m_p/m_s)*(R_p/a)^5'].iloc[0] * (-27 * np.pi /2)
    factor_25 = df[df['System'] == planet_25]['(m_p/m_s)*(R_p/a)^5'].iloc[0] * (-27 * np.pi /2)
    factor_26 = df[df['System'] == planet_26]['(m_p/m_s)*(R_p/a)^5'].iloc[0] * (-27 * np.pi /2)
    factor_27 = df[df['System'] == planet_27]['(m_p/m_s)*(R_p/a)^5'].iloc[0] * (-27 * np.pi /2)
    factor_28 = df[df['System'] == planet_28]['(m_p/m_s)*(R_p/a)^5'].iloc[0] * (-27 * np.pi /2)
    factor_29 = df[df['System'] == planet_29]['(m_p/m_s)*(R_p/a)^5'].iloc[0] * (-27 * np.pi /2)
    factor_30 = df[df['System'] == planet_30]['(m_p/m_s)*(R_p/a)^5'].iloc[0] * (-27 * np.pi /2)

    data_1 = pd.read_csv(path_1 + f'{planet_1}.csv')
    data_2 = pd.read_csv(path_2 + f'{planet_2}.csv')
    data_3 = pd.read_csv(path_3 + f'{planet_3}.csv')
    data_4 = pd.read_csv(path_4 + f'{planet_4}.csv')
    data_5 = pd.read_csv(path_5 + f'{planet_5}.csv')
    data_6 = pd.read_csv(path_6 + f'{planet_6}.csv')
    data_7 = pd.read_csv(path_7 + f'{planet_7}.csv')
    data_8 = pd.read_csv(path_8 + f'{planet_8}.csv')
    data_9 = pd.read_csv(path_9 + f'{planet_9}.csv')
    data_10 = pd.read_csv(path_10 + f'{planet_10}.csv')
    data_11 = pd.read_csv(path_11 + f'{planet_11}.csv')
    data_12 = pd.read_csv(path_12 + f'{planet_12}.csv')
    data_13 = pd.read_csv(path_13 + f'{planet_13}.csv')
    data_14 = pd.read_csv(path_14 + f'{planet_14}.csv')
    data_15 = pd.read_csv(path_15 + f'{planet_15}.csv')
    data_16 = pd.read_csv(path_16 + f'{planet_16}.csv')
    data_17 = pd.read_csv(path_17 + f'{planet_17}.csv')
    data_18 = pd.read_csv(path_18 + f'{planet_18}.csv')
    data_19 = pd.read_csv(path_19 + f'{planet_19}.csv')
    data_20 = pd.read_csv(path_20 + f'{planet_20}.csv')
    data_21 = pd.read_csv(path_21 + f'{planet_21}.csv')
    data_22 = pd.read_csv(path_22 + f'{planet_22}.csv')
    data_23 = pd.read_csv(path_23 + f'{planet_23}.csv')
    data_24 = pd.read_csv(path_24 + f'{planet_24}.csv')
    data_25 = pd.read_csv(path_25 + f'{planet_25}.csv')
    data_26 = pd.read_csv(path_26 + f'{planet_26}.csv')
    data_27 = pd.read_csv(path_27 + f'{planet_27}.csv')
    data_28 = pd.read_csv(path_28 + f'{planet_28}.csv')
    data_29 = pd.read_csv(path_29 + f'{planet_29}.csv')
    data_30 = pd.read_csv(path_30 + f'{planet_30}.csv')

    period_1 = data_1['Period'].iloc[0]
    t1 = data_1['Mid-point']
    uncertainty_1 = data_1['Uncertainty']

    period_2 = data_2['Period'].iloc[0]
    t2 = data_2['Mid-point']
    uncertainty_2 = data_2['Uncertainty']

    period_3 = data_3['Period'].iloc[0]
    t3 = data_3['Mid-point']
    uncertainty_3 = data_3['Uncertainty']

    period_4 = data_4['Period'].iloc[0]
    t4 = data_4['Mid-point']
    uncertainty_4 = data_4['Uncertainty']

    period_5 = data_5['Period'].iloc[0]
    t5 = data_5['Mid-point']
    uncertainty_5 = data_5['Uncertainty']

    period_6 = data_6['Period'].iloc[0]
    t6 = data_6['Mid-point']
    uncertainty_6 = data_6['Uncertainty']

    period_7 = data_7['Period'].iloc[0]
    t7 = data_7['Mid-point']
    uncertainty_7 = data_7['Uncertainty']

    period_8 = data_8['Period'].iloc[0]
    t8 = data_8['Mid-point']
    uncertainty_8 = data_8['Uncertainty']

    period_9 = data_9['Period'].iloc[0]
    t9 = data_9['Mid-point']
    uncertainty_9 = data_9['Uncertainty']

    period_10 = data_10['Period'].iloc[0]
    t10 = data_10['Mid-point']
    uncertainty_10 = data_10['Uncertainty']

    period_11 = data_11['Period'].iloc[0]
    t11 = data_11['Mid-point']
    uncertainty_11 = data_11['Uncertainty']

    period_12 = data_12['Period'].iloc[0]
    t12 = data_12['Mid-point']
    uncertainty_12 = data_12['Uncertainty']

    period_13 = data_13['Period'].iloc[0]
    t13 = data_13['Mid-point']
    uncertainty_13 = data_13['Uncertainty']

    period_14 = data_14['Period'].iloc[0]
    t14 = data_14['Mid-point']
    uncertainty_14 = data_14['Uncertainty']

    period_15 = data_15['Period'].iloc[0]
    t15 = data_15['Mid-point']
    uncertainty_15 = data_15['Uncertainty']

    period_16 = data_16['Period'].iloc[0]
    t16 = data_16['Mid-point']
    uncertainty_16 = data_16['Uncertainty']

    period_17 = data_17['Period'].iloc[0]
    t17 = data_17['Mid-point']
    uncertainty_17 = data_17['Uncertainty']

    period_18 = data_18['Period'].iloc[0]
    t18 = data_18['Mid-point']
    uncertainty_18 = data_18['Uncertainty']

    period_19 = data_19['Period'].iloc[0]
    t19 = data_19['Mid-point']
    uncertainty_19 = data_19['Uncertainty']

    period_20 = data_20['Period'].iloc[0]
    t20 = data_20['Mid-point']
    uncertainty_20 = data_20['Uncertainty']

    period_21 = data_21['Period'].iloc[0]
    t21 = data_21['Mid-point']
    uncertainty_21 = data_21['Uncertainty']

    period_22 = data_22['Period'].iloc[0]
    t22 = data_22['Mid-point']
    uncertainty_22 = data_22['Uncertainty']

    period_23 = data_23['Period'].iloc[0]
    t23 = data_23['Mid-point']
    uncertainty_23 = data_23['Uncertainty']

    period_24 = data_24['Period'].iloc[0]
    t24 = data_24['Mid-point']
    uncertainty_24 = data_24['Uncertainty']

    period_25 = data_25['Period'].iloc[0]
    t25 = data_25['Mid-point']
    uncertainty_25 = data_25['Uncertainty']

    period_26 = data_26['Period'].iloc[0]
    t26 = data_26['Mid-point']
    uncertainty_26 = data_26['Uncertainty']

    period_27 = data_27['Period'].iloc[0]
    t27 = data_27['Mid-point']
    uncertainty_27 = data_27['Uncertainty']

    period_28 = data_28['Period'].iloc[0]
    t28 = data_28['Mid-point']
    uncertainty_28 = data_28['Uncertainty']

    period_29 = data_29['Period'].iloc[0]
    t29 = data_29['Mid-point']
    uncertainty_29 = data_29['Uncertainty']

    period_30 = data_30['Period'].iloc[0]
    t30 = data_30['Mid-point']
    uncertainty_30 = data_30['Uncertainty']


    t_arr = t1.to_numpy()
    a0 = 0.5 * (np.min(t_arr) + np.max(t_arr))
    t0_1 = t_arr.flat[np.abs(t_arr - a0).argmin()]

    t_arr = t2.to_numpy()
    a0 = 0.5 * (np.min(t_arr) + np.max(t_arr))
    t0_2 = t_arr.flat[np.abs(t_arr - a0).argmin()]


    t_arr = t3.to_numpy()
    a0 = 0.5 * (np.min(t_arr) + np.max(t_arr))
    t0_3 = t_arr.flat[np.abs(t_arr - a0).argmin()]

    t_arr = t4.to_numpy()
    a0 = 0.5 * (np.min(t_arr) + np.max(t_arr))
    t0_4 = t_arr.flat[np.abs(t_arr - a0).argmin()]

    t_arr = t5.to_numpy()
    a0 = 0.5 * (np.min(t_arr) + np.max(t_arr))
    t0_5 = t_arr.flat[np.abs(t_arr - a0).argmin()]

    t_arr = t6.to_numpy()
    a0 = 0.5 * (np.min(t_arr) + np.max(t_arr))
    t0_6 = t_arr.flat[np.abs(t_arr - a0).argmin()]

    t_arr = t7.to_numpy()
    a0 = 0.5 * (np.min(t_arr) + np.max(t_arr))
    t0_7 = t_arr.flat[np.abs(t_arr - a0).argmin()]

    t_arr = t8.to_numpy()
    a0 = 0.5 * (np.min(t_arr) + np.max(t_arr))
    t0_8 = t_arr.flat[np.abs(t_arr - a0).argmin()]

    t_arr = t9.to_numpy()
    a0 = 0.5 * (np.min(t_arr) + np.max(t_arr))
    t0_9 = t_arr.flat[np.abs(t_arr - a0).argmin()]

    t_arr = t10.to_numpy()
    a0 = 0.5 * (np.min(t_arr) + np.max(t_arr))
    t0_10 = t_arr.flat[np.abs(t_arr - a0).argmin()]

    t_arr = t11.to_numpy()
    a0 = 0.5 * (np.min(t_arr) + np.max(t_arr))
    t0_11 = t_arr.flat[np.abs(t_arr - a0).argmin()]

    t_arr = t12.to_numpy()
    a0 = 0.5 * (np.min(t_arr) + np.max(t_arr))
    t0_12 = t_arr.flat[np.abs(t_arr - a0).argmin()]

    t_arr = t13.to_numpy()
    a0 = 0.5 * (np.min(t_arr) + np.max(t_arr))
    t0_13 = t_arr.flat[np.abs(t_arr - a0).argmin()]

    t_arr = t14.to_numpy()
    a0 = 0.5 * (np.min(t_arr) + np.max(t_arr))
    t0_14 = t_arr.flat[np.abs(t_arr - a0).argmin()]

    t_arr = t15.to_numpy()
    a0 = 0.5 * (np.min(t_arr) + np.max(t_arr))
    t0_15 = t_arr.flat[np.abs(t_arr - a0).argmin()]

    t_arr = t16.to_numpy()
    a0 = 0.5 * (np.min(t_arr) + np.max(t_arr))
    t0_16 = t_arr.flat[np.abs(t_arr - a0).argmin()]

    t_arr = t17.to_numpy()
    a0 = 0.5 * (np.min(t_arr) + np.max(t_arr))
    t0_17 = t_arr.flat[np.abs(t_arr - a0).argmin()]

    t_arr = t18.to_numpy()
    a0 = 0.5 * (np.min(t_arr) + np.max(t_arr))
    t0_18 = t_arr.flat[np.abs(t_arr - a0).argmin()]

    t_arr = t19.to_numpy()
    a0 = 0.5 * (np.min(t_arr) + np.max(t_arr))
    t0_19 = t_arr.flat[np.abs(t_arr - a0).argmin()]

    t_arr = t20.to_numpy()
    a0 = 0.5 * (np.min(t_arr) + np.max(t_arr))
    t0_20 = t_arr.flat[np.abs(t_arr - a0).argmin()]

    t_arr = t21.to_numpy()
    a0 = 0.5 * (np.min(t_arr) + np.max(t_arr))
    t0_21 = t_arr.flat[np.abs(t_arr - a0).argmin()]

    t_arr = t22.to_numpy()
    a0 = 0.5 * (np.min(t_arr) + np.max(t_arr))
    t0_22 = t_arr.flat[np.abs(t_arr - a0).argmin()]

    t_arr = t23.to_numpy()
    a0 = 0.5 * (np.min(t_arr) + np.max(t_arr))
    t0_23 = t_arr.flat[np.abs(t_arr - a0).argmin()]

    t_arr = t24.to_numpy()
    a0 = 0.5 * (np.min(t_arr) + np.max(t_arr))
    t0_24 = t_arr.flat[np.abs(t_arr - a0).argmin()]

    t_arr = t25.to_numpy()
    a0 = 0.5 * (np.min(t_arr) + np.max(t_arr))
    t0_25 = t_arr.flat[np.abs(t_arr - a0).argmin()]

    t_arr = t26.to_numpy()
    a0 = 0.5 * (np.min(t_arr) + np.max(t_arr))
    t0_26 = t_arr.flat[np.abs(t_arr - a0).argmin()]

    t_arr = t27.to_numpy()
    a0 = 0.5 * (np.min(t_arr) + np.max(t_arr))
    t0_27 = t_arr.flat[np.abs(t_arr - a0).argmin()]

    t_arr = t28.to_numpy()
    a0 = 0.5 * (np.min(t_arr) + np.max(t_arr))
    t0_28 = t_arr.flat[np.abs(t_arr - a0).argmin()]


    t_arr = t29.to_numpy()
    a0 = 0.5 * (np.min(t_arr) + np.max(t_arr))
    t0_29 = t_arr.flat[np.abs(t_arr - a0).argmin()]


    t_arr = t30.to_numpy()
    a0 = 0.5 * (np.min(t_arr) + np.max(t_arr))
    t0_30 = t_arr.flat[np.abs(t_arr - a0).argmin()]




    # orbit numbers for 30 planets
    orbit_number_1 = fold_time(t1, t0_1, period_1)
    orbit_number_2 = fold_time(t2, t0_2, period_2)
    orbit_number_3 = fold_time(t3, t0_3, period_3)
    orbit_number_4 = fold_time(t4, t0_4, period_4)
    orbit_number_5 = fold_time(t5, t0_5, period_5)
    orbit_number_6 = fold_time(t6, t0_6, period_6)
    orbit_number_7 = fold_time(t7, t0_7, period_7)
    orbit_number_8 = fold_time(t8, t0_8, period_8)
    orbit_number_9 = fold_time(t9, t0_9, period_9)
    orbit_number_10 = fold_time(t10, t0_10, period_10)
    orbit_number_11 = fold_time(t11, t0_11, period_11)
    orbit_number_12 = fold_time(t12, t0_12, period_12)
    orbit_number_13 = fold_time(t13, t0_13, period_13)
    orbit_number_14 = fold_time(t14, t0_14, period_14)
    orbit_number_15 = fold_time(t15, t0_15, period_15)
    orbit_number_16 = fold_time(t16, t0_16, period_16)
    orbit_number_17 = fold_time(t17, t0_17, period_17)
    orbit_number_18 = fold_time(t18, t0_18, period_18)
    orbit_number_19 = fold_time(t19, t0_19, period_19)
    orbit_number_20 = fold_time(t20, t0_20, period_20)
    orbit_number_21 = fold_time(t21, t0_21, period_21)
    orbit_number_22 = fold_time(t22, t0_22, period_22)
    orbit_number_23 = fold_time(t23, t0_23, period_23)
    orbit_number_24 = fold_time(t24, t0_24, period_24)
    orbit_number_25 = fold_time(t25, t0_25, period_25)
    orbit_number_26 = fold_time(t26, t0_26, period_26)
    orbit_number_27 = fold_time(t27, t0_27, period_27)
    orbit_number_28 = fold_time(t28, t0_28, period_28)
    orbit_number_29 = fold_time(t29, t0_29, period_29)
    orbit_number_30 = fold_time(t30, t0_30, period_30)





    q = 10000

    nll = lambda *args: -lnlike(*args) 
    initial = np.array([period_1, t0_1,
                        period_2, t0_2,
                        period_3, t0_3,
                        period_4, t0_4,
                        period_5, t0_5,
                        period_6, t0_6,
                        period_7, t0_7,
                        period_8, t0_8,
                        period_9, t0_9,
                        period_10, t0_10,
                        period_11, t0_11,
                        period_12, t0_12,
                        period_13, t0_13,
                        period_14, t0_14,
                        period_15, t0_15,
                        period_16, t0_16,
                        period_17, t0_17,
                        period_18, t0_18,
                        period_19, t0_19,
                        period_20, t0_20,
                        period_21, t0_21,
                        period_22, t0_22,
                        period_23, t0_23,
                        period_24, t0_24,
                        period_25, t0_25,
                        period_26, t0_26,
                        period_27, t0_27,
                        period_28, t0_28,
                        period_29, t0_29,
                        period_30, t0_30,
                        q]) 




    soln = minimize(nll, initial, args=(orbit_number_1, t1, uncertainty_1,
                                        orbit_number_2, t2, uncertainty_2,
                                        orbit_number_3, t3, uncertainty_3,
                                        orbit_number_4, t4, uncertainty_4,
                                        orbit_number_5, t5, uncertainty_5,
                                        orbit_number_6, t6, uncertainty_6,
                                        orbit_number_7, t7, uncertainty_7,
                                        orbit_number_8, t8, uncertainty_8,
                                        orbit_number_9, t9, uncertainty_9,
                                        orbit_number_10, t10, uncertainty_10,
                                        orbit_number_11, t11, uncertainty_11,
                                        orbit_number_12, t12, uncertainty_12,
                                        orbit_number_13, t13, uncertainty_13,
                                        orbit_number_14, t14, uncertainty_14,
                                        orbit_number_15, t15, uncertainty_15,
                                        orbit_number_16, t16, uncertainty_16,
                                        orbit_number_17, t17, uncertainty_17,
                                        orbit_number_18, t18, uncertainty_18,
                                        orbit_number_19, t19, uncertainty_19,
                                        orbit_number_20, t20, uncertainty_20,
                                        orbit_number_21, t21, uncertainty_21,
                                        orbit_number_22, t22, uncertainty_22,
                                        orbit_number_23, t23, uncertainty_23,
                                        orbit_number_24, t24, uncertainty_24,
                                        orbit_number_25, t25, uncertainty_25,
                                        orbit_number_26, t26, uncertainty_26,
                                        orbit_number_27, t27, uncertainty_27,
                                        orbit_number_28, t28, uncertainty_28,
                                        orbit_number_29, t29, uncertainty_29,
                                        orbit_number_30, t30, uncertainty_30,
                                        factor_1,
                                        factor_2,
                                        factor_3,
                                        factor_4,
                                        factor_5,
                                        factor_6,
                                        factor_7,
                                        factor_8,
                                        factor_9,
                                        factor_10,
                                        factor_11,
                                        factor_12,
                                        factor_13,
                                        factor_14,
                                        factor_15,
                                        factor_16,
                                        factor_17,
                                        factor_18,
                                        factor_19,
                                        factor_20,
                                        factor_21,
                                        factor_22,
                                        factor_23,
                                        factor_24,
                                        factor_25,
                                        factor_26,
                                        factor_27,
                                        factor_28,
                                        factor_29,
                                        factor_30
        ), method='Nelder-Mead')  





    print(soln.x) 






main()










