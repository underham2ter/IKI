import os
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from misc_functions import *
import disort

def call_disort(lines_path, input_path, profile_filename, aerosol_filename, gas_ext_path, gas_ext_file, hmin=0, emis_ang=0, sf=1, emis=0.95, ic=0, dh=0):

    Loshmidt = 2.69e19

    GAS_EXT_PATH = gas_ext_path
    GAS_EXT_FILE = gas_ext_file
    
    PROFILE_FILENAME = profile_filename
    AEROSOL_FILENAME = aerosol_filename
    
    LINES_PATH = lines_path
    INPUT_PATH = input_path

    CORE_NAME = GAS_EXT_FILE[8:-4]
    SF = str(sf).replace('.', '_')
    IC = str(ic).replace('.', '_')
    
    OUTPUT_PATH = f'disort_output/{CORE_NAME}/'
    OUTPUT_FILENAME = f'disort_result_emis_{emis}_ang_{emis_ang}_hmin_{hmin}_sf_{SF}_ic_{IC}_{CORE_NAME}.csv'
    os.makedirs(OUTPUT_PATH, exist_ok=True)

    pressure, height, concentration, temperature = np.loadtxt(INPUT_PATH + PROFILE_FILENAME, skiprows = 1, unpack=True)
    
    height_corrected = np.arange(hmin, 101, 1)
    pressure_corrected = np.exp(interp1d(height, np.log(pressure), kind='linear', fill_value='extrapolate')(height_corrected))
    concentration_corrected = np.exp(interp1d(height, np.log(concentration), kind='linear', fill_value='extrapolate')(height_corrected))
    temperature_corrected = interp1d(height, temperature, kind='linear', fill_value='extrapolate')(height_corrected)
    
    atm = np.column_stack([
        pressure_corrected, 
        height_corrected, 
        concentration_corrected, 
        temperature_corrected
    ])
    
    gas_abs_data = pd.read_csv(GAS_EXT_PATH + GAS_EXT_FILE, header=None)
    gas_abs_data = np.array(gas_abs_data)
    
    gas_ext = gas_abs_data[1:]
    wavenum = gas_abs_data[0]
    
    # gas_ext = interp1d(height, gas_ext, kind='linear', fill_value='extrapolate', axis=0)(height_corrected)
    gas_ext = gas_ext[hmin+3:]
    gas_ext = (gas_ext.T + ic * 1e-9/Loshmidt**2 * atm[:, 2]**2).T
    
    # Оптические свойства аэрозолей
    data = np.loadtxt(INPUT_PATH + AEROSOL_FILENAME)
    
    # Расчет рассеяния и поглощения аэрозолей
    modeH = [1, 21, 22, 3]
    aer_wavenum = data[::4, 0]
    wavelen = data[::4, 1]
    
    modes = data[:, 2].reshape(len(aer_wavenum), 4).T
    Cext = data[:, 3].reshape(len(aer_wavenum), 4)
    Csca = data[:, 6].reshape(len(aer_wavenum), 4)
    
    Galb = data[:, 4].reshape(len(aer_wavenum), 4)
    Gass = data[:, 5].reshape(len(aer_wavenum), 4)
    
    
    nphi = 51
    nL = 101
    MILC = data[:, 7 + nphi * 2:].reshape(len(aer_wavenum), 4, nL)
    mileg = np.zeros((4, nL, len(aer_wavenum)))
    
    
    for i in range(4):
        mileg[:, :, i] = MILC[i, :, :]
    
    zb = np.array([49, 65, 49, 49])
    zb[:2] += dh
    zc = np.array([16, 1, 11, 8])
    N0 = np.array([193.5, 100, 50, 14])
    N0[1:] *= sf
    Hup = np.array([3.5, 3.5, 1.0, 1.0])
    Hlo = np.array([1.0, 3.0, 0.1, 0.5])
    
    N = np.zeros((len(atm[:, 1]), 4))
    
    for i in range(4):
        idx1 = atm[:, 1] > (zb[i] + zc[i])
        idx2 = (atm[:, 1] <= (zb[i] + zc[i])) & (atm[:, 1] >= zb[i])
        idx3 = atm[:, 1] < zb[i]
        
        N[idx1, i] = N0[i] * np.exp(-(atm[idx1, 1] - (zb[i] + zc[i])) / Hup[i])
        N[idx2, i] = N0[i]
        N[idx3, i] = N0[i] * np.exp(-(zb[i] - atm[idx3, 1]) / Hlo[i])
    
    aer_ext = np.dot(Cext, N.T) / 1000  # Cext [micron^2] * N [1/cm^3] * 10^(-3) --> km^(-1)
    aer_sca = np.dot(Csca, N.T) / 1000
    
    aer_ext = aer_ext.T
    aer_sca = aer_sca.T
    
    gfmat = Gass
    
    scamat = np.zeros((len(wavelen), len(atm[:, 1]), 4))
    for i in range(4):
        scamat[:, :, i] = np.dot(Csca[:, i].reshape(-1, 1), N[:, i].reshape(1, -1)) / 1000
    
    aer_ssalb = aer_sca / aer_ext
    aer_ssalb[aer_ext <= 0] = 0
    
    # Вычисление длины волны
    lambda_vals = 1e4 / wavenum
    
    # Инициализация массива для рассеяния Рэлея
    rayleigh = np.zeros((len(atm[:, 2]), len(wavenum)))
    
    # Цикл для вычисления рассеяния Рэлея для каждого значения в атмосфере
    for i in range(len(atm[:, 1])):
        mus = (6.99100e-2 / (166.175000 - lambda_vals**(-2)) +
               1.44720e-3 / (79.609000 - lambda_vals**(-2)) +
               6.42941e-5 / (56.306400 - lambda_vals**(-2)) +
               5.21306e-5 / (46.019600 - lambda_vals**(-2)) +
               1.46847e-6 / (0.0584738 - lambda_vals**(-2)))
        
        n = 1 + mus * atm[i, 0] / atm[i, 3] * 273.15 / 1.01325
        
        rayleigh[i, :] = (24 * np.pi**3 / (lambda_vals * 1e-4)**4 / atm[i, 2]**2 *
                          ((n**2 - 1) / (n**2 + 2))**2 *
                          (1.1364 + 25.3e-4 / lambda_vals**2) * atm[i, 2])
    
    # Преобразование единиц
    rayleigh *= 1e5
    
    
    SFCTEMP = atm[0, 3]
    
    ZZ = atm[::-1, 1]
    TT = atm[::-1, 3]
    delta_zz = 1 #np.abs(np.diff(ZZ))
    
    GF = np.zeros_like(gas_ext)  # текущий геометрический фактор уровня
    
    aer_ext_interp = interp1d(aer_wavenum, aer_ext, axis=1, kind='linear', fill_value='extrapolate')(wavenum)
    aer_sca_interp = interp1d(aer_wavenum, aer_sca, axis=1, kind='linear', fill_value='extrapolate')(wavenum)
    
    TAUCL = gas_ext + aer_ext_interp + rayleigh  # текущий уровень экстинкции
    SCACL = aer_sca_interp + rayleigh            # текущий уровень рассеяния
    TAUCL = TAUCL[::-1, :]                       # однократное рассеяние (однородное поглощение)
    SCACL = SCACL[::-1, :]
    SSALB = SCACL / TAUCL
    
    for i in range(len(ZZ)-1, -1, -1):
        i1 = len(ZZ) - i - 1
    
        GF[i1,:] = (interp1d(aer_wavenum, gfmat[:, 0] * scamat[:, i, 0], kind='linear', fill_value='extrapolate')(wavenum) +
                     interp1d(aer_wavenum, gfmat[:, 1] * scamat[:, i, 1], kind='linear', fill_value='extrapolate')(wavenum) +
                     interp1d(aer_wavenum, gfmat[:, 2] * scamat[:, i, 2], kind='linear', fill_value='extrapolate')(wavenum) +
                     interp1d(aer_wavenum, gfmat[:, 3] * scamat[:, i, 3], kind='linear', fill_value='extrapolate')(wavenum))  / SCACL[i1,:]
    
    TAUSUM = (TAUCL[:-1, :] + TAUCL[1:, :]) / 2 #* np.repeat(delta_zz[np.newaxis, :], len(aer_wavenum), axis=0).T
    SCASUM = (SCACL[:-1, :] + SCACL[1:, :]) / 2 #* np.repeat(delta_zz[np.newaxis, :], len(aer_wavenum), axis=0).T
    
    SSALB = (SSALB[:-1, :] + SSALB[1:, :]) / 2
    GF = (GF[:-1, :] + GF[1:, :]) / 2
    
    NLYR = len(ZZ)-1
    
    UMU = [np.cos(np.radians(emis_ang))]
    
    # UMU = [-1.0, -0.9848, -0.9397, -0.8660, -0.7660, -0.6428, -0.5, 
    #        -0.3420, -0.1736, -0.0872, -.01, .01,  0.0872, 0.1736,
    #        0.3420,  0.5,  0.6428,  0.7660,  0.8660,  0.9397,  0.9848, 1.0]
    
    # UMU = np.array([-1.,-0.5,0.5,1.])
    
    PHI0 = 0
    PHI = [0]
    
    UMU0 = 1
    FISOT = 0.0
    
    IBCND = 0
    FBEAM = 0
    
    SFCEMIS = emis # эмиссионная способность поверхности
    LAMBER = True
    ALBEDO = 1. - SFCEMIS
    PLANK = True
    TTEMP = 0.
    TEMIS = 1.0
    
    ONLYFL = False
    USRTAU = False
    USRANG = True
    PRNT = np.array([False, False, False, False, False])
    # PRNT = np.array([True, True, True, True, True])
    
    NSTR = 16 # 32
    MAXMOM = 299
    MAXCLY = NLYR
    MAXUMU = len(UMU)
    MAXPHI = len(PHI)
    MAXULV = NLYR + 1
    
    ACCUR = 0
    
    IPHAS = [3]*NLYR
    
    UTAU = np.zeros(NLYR+1)
    
    EARTH_RADIUS = 6052.1
    DO_PSEUDO_SPHERE = True
    DELTAMPLUS = True
    
    H_LYR = ZZ
    
    uu_ = []
    
    for i, w_iter in enumerate(wavenum):
        [rfldir, rfldn, flup, dfdt, uavg, uu, albmed, trnmed, rhoq, rhou, rho_a, bemst, emust] =\
                                    disort.run(dtauc=TAUSUM[:,i], ssalb=SSALB[:,i], iphas=IPHAS, gg=GF[:,i],
                                                umu0=UMU0, phi0=PHI0, albedo=ALBEDO, fbeam=FBEAM,
                                                usrtau=USRTAU, utau=UTAU, usrang=USRANG, umu=UMU,
                                                phi=PHI, nstr=NSTR, maxmom=MAXMOM, lamber=LAMBER,
                                                onlyfl=ONLYFL, accur=ACCUR, plank=PLANK,
                                                temper=TT, wvnmlo=w_iter, wvnmhi=w_iter,
                                                ibcnd=IBCND, fisot=FISOT, btemp=SFCTEMP, ttemp=TTEMP,
                                                temis=TEMIS, prnt=PRNT, earth_radius=EARTH_RADIUS, 
                                                h_lyr=H_LYR, deltamplus=DELTAMPLUS, 
                                                do_pseudo_sphere=DO_PSEUDO_SPHERE)
        uu_.append(uu)
    
    uu_ = np.array(uu_)

    return wavenum, uu_[:, 0, 0, 0]

    # output = np.vstack((wavenum, uu_[:, :, 0, 0].T)) 
    
    #          [wavenum]
    # UMU[0]   [intensity]
    # UMU[1]   [intensity]
    # ...
    # UMU[-1]  [intensity]
    
    # np.savetxt(OUTPUT_PATH + OUTPUT_FILENAME, output, delimiter=',')

def interpolate_spectrum(sf, emis_ang, topo, detector_num, interp_func, spicav_wavelength, disort_wavenum):
    sf = sf[0] if isinstance(sf, (list, np.ndarray)) else sf

    STEP = disort_wavenum[1] - disort_wavenum[0]

    interp_disort_spectrum = interp_func([sf, topo, emis_ang])[0]

    conv_disort_spectrum = conv_spicavir_matrix_sw_lw(interp_disort_spectrum, 
                                                     STEP, 
                                                     CONV_PATH, 
                                                     'SW', 
                                                     1, 
                                                     detector_num)

    conv_spectrum_wavelength_lattice = conv_disort_spectrum * disort_wavenum**2 / 1e4
    disort_wavelength = 1e4 / disort_wavenum
    
    final_lambda_lattice, final_disort_spectrum, spicav_mask =\
        integrate_model_on_data_lattice(spicav_wavelength,
                                        disort_wavelength,
                                        conv_spectrum_wavelength_lattice)
    
    return final_lambda_lattice, final_disort_spectrum, spicav_mask

# Функция потерь для оптимизации (например, среднеквадратичная ошибка)
def loss_function(sf, emis_ang, topo, detector_num, interp_func, spicav_spectrum, spicav_wavelength, disort_wavenum):

    final_lambda_lattice, final_disort_spectrum, spicav_mask =\
        interpolate_spectrum(sf, emis_ang, topo, detector_num, interp_func, spicav_wavelength, disort_wavenum)

    opt_cutoff = int(2*len(final_disort_spectrum)/3)

    loss_metric = np.mean((final_disort_spectrum - spicav_spectrum[::-1][spicav_mask]) ** 2)
    # loss_metric = np.mean((final_disort_spectrum[opt_cutoff:] - spicav_spectrum[::-1][spicav_mask][opt_cutoff:]) ** 2)
    
    return loss_metric