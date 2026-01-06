import numpy as np
import mpmath as mpm
from .RELAP_v4 import *




# single porosity

def simulateSinglePorosity(
        mean_residence_time,
        peclet_number,
        frac_retard,
        time_points,
        bckgrnd_conc,
        inj_concs,
        inj_durs,
        recRatio=0,
        wsCoef=0,
        dps = 8
):
    
    mpm.mp.dps = dps
    gw_inf = GroundWaterInfiniteMatrixSolution(
        mean_residence_time,
        peclet_number,
        0,
        0,
        fracture_retardation=frac_retard,
        matrix_retardation=1
    )

    tracer_injection = Input_Pulses_of_Tracer(
        background_concentration=bckgrnd_conc,
        injection_concentrations=inj_concs, # mg/L
        injection_durations=inj_durs

    )

    reciculation = Recirculation(recRatio)
    if wsCoef > 0:
        wellbore_storage = WellboreStorage(wsCoef)
    else:
        wellbore_storage = lambda x: 1
    relap_instance = RELAP_Modifed(
        Ground_Water=gw_inf,
        Input_Instance=tracer_injection,
        recirculation=reciculation,
        wellbore_storage_node=wellbore_storage
    )

    conc_values = Simulate_RELAP_Relative(
        relap_instance,
        time_points=time_points,
        background_concentration=bckgrnd_conc
    )

    return np.array(conc_values, dtype=np.float64)

def simulateDualPorosity(
        mean_residence_time,
        peclet_number,
        frac_retard,
        time_points,
        bckgrnd_conc,
        inj_concs,
        inj_durs,
        recRatio,
        dualPorosity_param,
        matrix_retardation,
        wsCoef=0,
        delay_time=0,
        dps = 8
):
    
    mpm.mp.dps = dps
    gw_inf = GroundWaterInfiniteMatrixSolution(
        mean_residence_time,
        peclet_number,
        dualPorosity_param,
        0,
        fracture_retardation=frac_retard,
        matrix_retardation=matrix_retardation
    )

    tracer_injection = Input_Pulses_of_Tracer(
        background_concentration=bckgrnd_conc,
        injection_concentrations=inj_concs, # mg/L
        injection_durations=inj_durs

    )

    reciculation = Recirculation(recRatio)
    pipelinedelay = PipelineDelay(delay_time=delay_time)
    if wsCoef > 0:
        wellbore_storage = WellboreStorage(wsCoef)
    else:
        wellbore_storage = lambda x: 1

    relap_instance = RELAP_Modifed(
        Ground_Water=gw_inf,
        Input_Instance=tracer_injection,
        recirculation=reciculation,
        pipeline_delay_node=pipelinedelay,
        wellbore_storage_node=wellbore_storage
    )

    conc_values = Simulate_RELAP_Relative(
        relap_instance,
        time_points=time_points,
        background_concentration=bckgrnd_conc
    )

    return np.array(conc_values, dtype=np.float64)

def simulateDualPorosityInf(
        mean_residence_time,
        peclet_number,
        matr_diff,
        poro_ratio,
        frac_retard,
        mtrx_retard,
        time_points,
        bckgrnd_conc,
        inj_concs,
        inj_durs,
        dps = 8
):
    
    mpm.mp.dps = dps
    gw_inf = GroundWaterInfiniteMatrixSolution(
        mean_residence_time,
        peclet_number,
        matr_diff,
        poro_ratio,
        0,
        fracture_retardation=frac_retard,
        matrix_retardation=mtrx_retard
    )

    tracer_injection = Input_Pulses_of_Tracer(
        background_concentration=bckgrnd_conc,
        injection_concentrations=inj_concs, # mg/L
        injection_durations=inj_durs

    )

    relap_instance = RELAP_Modifed(
        Ground_Water=gw_inf,
        Input_Isntance=tracer_injection
    )

    conc_values = Simulate_RELAP_Relative(
        relap_instance,
        time_points=time_points,
        background_concentration=bckgrnd_conc
    )

    return np.array(conc_values, dtype=np.float64)


def simulateDualPorosityfinite(
        mean_residence_time,
        peclet_number,
        matr_diff,
        poro_ratio,
        wtr_ratio,
        frac_retard,
        mtrx_retard,
        time_points,
        bckgrnd_conc,
        inj_concs,
        inj_durs,
        dps = 8
):
    
    mpm.mp.dps = dps
    gw_inf = GroundWaterFiniteMatrixSolution_2(
        mean_residence_time,
        peclet_number,
        matr_diff,
        poro_ratio,
        wtr_ratio,
        0,
        fracture_retardation=frac_retard,
        matrix_retardation=mtrx_retard
    )

    tracer_injection = Input_Pulses_of_Tracer(
        background_concentration=bckgrnd_conc,
        injection_concentrations=inj_concs, # mg/L
        injection_durations=inj_durs

    )

    relap_instance = RELAP_Modifed(
        Ground_Water=gw_inf,
        Input_Isntance=tracer_injection
    )

    conc_values = Simulate_RELAP_Relative(
        relap_instance,
        time_points=time_points,
        background_concentration=bckgrnd_conc
    )

    return np.array(conc_values, dtype=np.float64)


def simulateRoseNDSSinglePoro(
        mean_residence_time,
        peclet_number,
        dps=8
):
    
    frac_retard = 1.0
    bckgrnd_conc = 0.0
    inj_concs = np.array([7.0, 0], dtype=float)
    inj_durs = np.cumsum(np.array([1.5]))
    time_points = np.array(
        [  0.000001    ,   1.12625698,   2.41340782,   3.5396648 ,
           4.97640704,   5.91694907,   7.16864121,   8.41146867,
           9.20175447,  10.63871832,  11.44606861,  13.39718362,
          15.34807701,  17.46097234,  19.40144974,  21.50348585,
          22.79484741,  23.92398541,  26.02247565,  27.95763424,
          30.05302184,  31.98685073,  32.95354356,  35.20450621,
          36.49564615,  40.52043031,  45.51257214,  47.28417748,
          48.89422411,  50.98695231,  53.07746433,  55.00929867,
          56.94113302,  59.1938686 ,  61.44638257,  65.4696154 ,
          70.62043492,  73.03672376,  74.96745002,  77.21996399,
          79.15069024,  81.0814165 ,  83.17303661,  85.42555058,
          89.4485618 ,  96.52811302,  98.94152085, 101.03336257,
         102.96475368, 106.98732167, 109.40095111, 113.58441295,
         119.37570525, 121.30643151, 123.237601  , 125.0074334 ,
         127.09883189, 128.87065885], dtype=float)

    
    concs = simulateSinglePorosity(
        mean_residence_time,
        peclet_number,
        frac_retard,
        time_points,
        bckgrnd_conc,
        inj_concs,
        inj_durs,
        dps = dps
)
    
    return concs


def simulateRoseNDSDualPoroInf(
        mean_residence_time,
        peclet_number,
        matr_diff,
        poro_ratio,
        dps=8
):
    
    frac_retard = 1.0
    mtrx_retard = 1.0
    bckgrnd_conc = 0.0
    inj_concs = np.array([7.0, 0], dtype=float)
    inj_durs = np.cumsum(np.array([1.5]))
    time_points = np.array(
        [  0.000001    ,   1.12625698,   2.41340782,   3.5396648 ,
           4.97640704,   5.91694907,   7.16864121,   8.41146867,
           9.20175447,  10.63871832,  11.44606861,  13.39718362,
          15.34807701,  17.46097234,  19.40144974,  21.50348585,
          22.79484741,  23.92398541,  26.02247565,  27.95763424,
          30.05302184,  31.98685073,  32.95354356,  35.20450621,
          36.49564615,  40.52043031,  45.51257214,  47.28417748,
          48.89422411,  50.98695231,  53.07746433,  55.00929867,
          56.94113302,  59.1938686 ,  61.44638257,  65.4696154 ,
          70.62043492,  73.03672376,  74.96745002,  77.21996399,
          79.15069024,  81.0814165 ,  83.17303661,  85.42555058,
          89.4485618 ,  96.52811302,  98.94152085, 101.03336257,
         102.96475368, 106.98732167, 109.40095111, 113.58441295,
         119.37570525, 121.30643151, 123.237601  , 125.0074334 ,
         127.09883189, 128.87065885], dtype=float)

    
    concs = simulateDualPorosityInf(
        mean_residence_time,
        peclet_number,
        matr_diff,
        poro_ratio,
        frac_retard,
        mtrx_retard,
        time_points,
        bckgrnd_conc,
        inj_concs,
        inj_durs,
        dps=dps
        
)
    
    return concs

def simulateRoseNDSDualPoroFinite(
        mean_residence_time,
        peclet_number,
        matr_diff,
        poro_ratio,
        wtr_ratio,
        dps=8
):
    
    frac_retard = 1.0
    mtrx_retard = 1.0
    bckgrnd_conc = 0.0
    inj_concs = np.array([7.0, 0], dtype=float)
    inj_durs = np.cumsum(np.array([1.5]))
    time_points = np.array(
        [  0.000001    ,   1.12625698,   2.41340782,   3.5396648 ,
           4.97640704,   5.91694907,   7.16864121,   8.41146867,
           9.20175447,  10.63871832,  11.44606861,  13.39718362,
          15.34807701,  17.46097234,  19.40144974,  21.50348585,
          22.79484741,  23.92398541,  26.02247565,  27.95763424,
          30.05302184,  31.98685073,  32.95354356,  35.20450621,
          36.49564615,  40.52043031,  45.51257214,  47.28417748,
          48.89422411,  50.98695231,  53.07746433,  55.00929867,
          56.94113302,  59.1938686 ,  61.44638257,  65.4696154 ,
          70.62043492,  73.03672376,  74.96745002,  77.21996399,
          79.15069024,  81.0814165 ,  83.17303661,  85.42555058,
          89.4485618 ,  96.52811302,  98.94152085, 101.03336257,
         102.96475368, 106.98732167, 109.40095111, 113.58441295,
         119.37570525, 121.30643151, 123.237601  , 125.0074334 ,
         127.09883189, 128.87065885], dtype=float)

    
    concs = simulateDualPorosityfinite(
        mean_residence_time,
        peclet_number,
        matr_diff,
        poro_ratio,
        wtr_ratio,
        frac_retard,
        mtrx_retard,
        time_points,
        bckgrnd_conc,
        inj_concs,
        inj_durs,
        dps = dps
)
    
    return concs