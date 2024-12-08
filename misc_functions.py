import numpy as np
from scipy.signal import convolve
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

# Глобальные параметры
c = 2.99792458e10 # cms^−1
PI = 3.1415926535898  # Число Пи
PLANCK = 6.62606876E-27  # Постоянная Планка, эрг*с
BOLTZ = 1.3806503E-16  # Постоянная Больцмана, эрг/К
CLIGHT = 2.99792458E+10  # Скорость света в вакууме, см/с
AVOGAD = 6.02214199E+23  # Число Авогадро, 1/моль
ALOSMT = 2.6867775E+19  # Постоянная Лошмидта, молекулы/см³ (стандартное условие)
GASCON = 8.314472E+07  # Универсальная газовая постоянная, эрг/(моль*К)
RADCN1 = 1.191042722E-12  # Первая радиационная постоянная, эрг*см²/с
RADCN2 = 1.4387752  # Вторая радиационная постоянная, см*К

def doppler_fwhm(nu, T, molm):

    return nu/c*np.sqrt(2*AVOGAD*BOLTZ*T*np.log(2)/molm)


def lorenz_fwhm(p, T, T_ref, gamma_self, gamma_air, p_self, n):

    return ((T_ref/T)**n)*(gamma_self*p_self + gamma_air*(p-p_self))
    

def compute_profile(x, y):
    t = np.complex128(y - 1j * x)
    u = t * t
    s = np.abs(x) + y
    
    p = np.zeros(len(x))

    # region 1 (s>=15)
    ff1 = np.where(s >= 15)
    # print(ff1)
    # print(len(ff1[0]))
    if len(ff1[0]) > 0:
        w4 = t[ff1] * 0.5641896 / (0.5 + u[ff1])
        p[ff1] = np.real(w4)

    # region 2 (5.5 <= s < 15)
    ff2 = np.where((s < 15) & (s >= 5.5))
    # print(ff2)
    if len(ff2[0]) > 0:
        w4 = t[ff2] * (1.410474 + u[ff2] * 0.5641896) / (0.75 + u[ff2] * (3.0 + u[ff2]))
        p[ff2] = np.real(w4)

    # region 3 (s < 5.5, y >= 0.195 * abs(x) - 0.176)
    ff3 = np.where((s < 5.5) & (y >= 0.195 * np.abs(x) - 0.176))
    # print(ff3)
    if len(ff3[0]) > 0:
        w4 = (16.4955 + t[ff3] * 20.20933 + u[ff3] * 11.96482 + t[ff3] * u[ff3] * 3.778987 +
              u[ff3] * u[ff3] * 0.5642236) / (16.4955 + t[ff3] * 38.82363 + u[ff3] * 39.27121 +
              t[ff3] * u[ff3] * 21.69274 + u[ff3] * u[ff3] * 6.699398 + t[ff3] * u[ff3] * u[ff3])
        p[ff3] = np.real(w4)

    # region 4 (s < 5.5, y < 0.195 * abs(x) - 0.176)
    ff4 = np.where((s < 5.5) & (y < 0.195 * np.abs(x) - 0.176))
    # print(ff4)
    if len(ff4[0]) > 0:
        w4 = (np.exp(u[ff4]) - t[ff4] * (36183.31 - u[ff4] * (3321.9905 - u[ff4] * 
                (1540.787 - u[ff4] * (219.0313 - u[ff4] * (35.76683 - u[ff4] * 
                (1.320522 - u[ff4] * 0.56419))))))) / (32066.6 - u[ff4] * 
                (24322.84 - u[ff4] * (9022.228 - u[ff4] * (2186.181 - u[ff4] * 
                (364.2191 - u[ff4] * (61.57037 - u[ff4] * (1.841439 - u[ff4])))))))
        p[ff4] = np.real(w4)

    return p 

def conv_spicavir_matrix_v2(spectrum, step, path, *args):
    """
    1D сверка, если сетка по числам волн равномерна, т.е. step = const
    spectrum -- спектр для свертки
    step       -- шаг сетки спектра (сетка возрастает!)
    """

    if len(args) < 2:
        file = 'psf_lw_all_O2_1270_desc.txt'
        data = np.loadtxt(path + file, skiprows=1)
        delta_nu = data[:, 0]
        profile = data[:, 2]
    else:
        delta_nu = args[0]
        profile = args[1]
    
    interp_func = interp1d(delta_nu, profile, kind='linear', fill_value=0, bounds_error=False)
    profile_interp = interp_func(np.arange(-100, 100 + step, step))
    
    norm = convolve(np.ones(len(spectrum)), profile_interp, mode='same')
    sp = convolve(spectrum, profile_interp, mode='same') / norm

    return sp

def conv_spicavir_matrix_sw_lw(spectrum, step, path, *args):
    """
    1D свертка, если сетка по числам волн равномерна, т.е. step = const
    spectrum -- спектр для свертки
    step       -- шаг сетки спектра (сетка возрастает!)
    
    args[0] -- канал SPICAV:
       'LW'  -- канал длинноволновый (может быть пропущен)
       'SW'  -- канал коротковолновый (!!!)
    args[1] -- ТОЛЬКО для SW канала:
       0.6   -- в микронах, эталонная длина волны: PSF зависит от частоты AOTF,
                то есть отличается для различных длин волн
    args[2] -- ТОЛЬКО для SW канала:
       0 или 1 -- номер детектора
    """

    if len(args) < 1:
        file = path + 'psf_lw_all_O2_1270_desc.txt'
        data = np.loadtxt(file)
        wdata = data[:, 0]
        fdata = data[:, 2]
    else:
        if args[0] == 'LW':
            file = path + 'psf_lw_all_O2_1270_desc.txt'
            data = np.loadtxt(file)
            wdata = data[:, 0]
            fdata = data[:, 2]
        elif args[0] == 'SW':
            file = path + 'psf_sw_all_editted.txt'
            data = np.loadtxt(file)
            dref = data[0, 1:]
            wref = data[1, 1:] / 1000  # из нм в микроны
            a = np.array([-4.9405101e-08, -5.0454785e-08])
            b = np.array([7.6969006e-02, 7.7358519e-02])
            c = np.array([-2.9822051e+02, -3.3244465e+02])
            
            # Расчет корней уравнения для выбранного детектора
            det_index = int(args[2])  # args[2] = 0 или 1
            r = np.roots([a[det_index], b[det_index], c[det_index] - 10**4 / args[1]])
            r = r[(r > 0.5e5) & (r < 3e5)]
            
            deltaf = data[2:, 0]
            wdata = (2 * a[det_index] * r + b[det_index]) * deltaf
            
            # Поиск индекса минимальной разницы между эталонной длиной волны и args[1]
            iref = np.argmin(np.abs(wref[::2] - args[1])) * 2
            iref = 1 + np.where((wref == wref[iref]) & (dref == args[2]))[0][0]
            
            fdata = data[2:, iref]
    
    # Интерполяция fdata для новых значений wdata
    interp_func = interp1d(wdata, fdata, kind='linear', fill_value=0, bounds_error=False)
    fdata_interp = interp_func(np.arange(-65, 65 + step, step))
    
    # Проведение свертки с полученным интерполированным fdata
    relsp = np.ones_like(spectrum)
    sp = convolve(spectrum, fdata_interp, mode='same') / convolve(relsp, fdata_interp, mode='same')

    return sp


def plot_for_selected_height(intensity, wavenum, angle_values):
    plt.figure()
    
    # Строим линии для каждого угла при выбранной высоте
    num_angles = intensity.shape[0]
    
    for angle_idx in range(12, num_angles, 3):
        intensities = intensity[angle_idx]
        label = r'$\cos \theta$' + f': {angle_values[angle_idx]}'
        plt.plot(wavenum, intensities*wavenum**2/1e4, label=label, linewidth=0.5)
        
    plt.xlabel(r'$\tilde{\nu}$, $cm^{-1}$')
    plt.ylabel(r'Intensity, $W/m^{2}/\mu m / sr$')
    plt.title(f'Intensity at the top boundary')
    plt.legend()
    plt.grid(True)
    plt.show()

def plot_for_grid(intensity, wavenum, angle_values):
    plt.figure()
    
    # Строим линии для каждого угла при выбранной высоте
    num_angles = intensity.shape[0]
    
    for angle_idx in range(0, num_angles, 1):
        intensities = intensity[angle_idx]
        label = r'param' + f': {angle_values[angle_idx]}'
        plt.plot(wavenum, intensities*wavenum**2/1e4, label=label, linewidth=0.5)
        
    plt.xlabel(r'$\tilde{\nu}$, $cm^{-1}$')
    plt.ylabel(r'Intensity')
    plt.title(f'Intensity at the top boundary')
    plt.legend()
    plt.grid(True)
    plt.show()


# def find_nearest(array,value):
#     idx = np.searchsorted(array, value, side="left")
#     if idx > 0 and (idx == len(array) or abs(value - array[idx-1]) < abs(value - array[idx])):
#         return idx-1
#     else:
#         return idx

def find_nearest(array, value):
    idx = np.searchsorted(-array, -value, side='left')
    if idx > 0 and (idx == len(array) or abs(value - array[idx-1]) < abs(value - array[idx])):
        return idx-1
    else:
        return idx

def integrate_model_on_data_lattice(data_lattice, model_lattice, spectra):
    spectra = np.array(spectra)
    
    if not np.all(np.diff(data_lattice) > 0):
        data_lattice = np.sort(data_lattice)
    
    if not np.all(np.diff(model_lattice) > 0):
        sorted_indices = np.argsort(model_lattice)
        model_lattice = model_lattice[sorted_indices]
        spectra = spectra[sorted_indices]
    
    # 1. Найдем средние точки между элементами массива data_lattice
    w_data_midpoints = (data_lattice[:-1] + data_lattice[1:]) / 2
    # 2. Создаем списки для результирующих частот и интегрированных значений
    w_data_result = []
    integrated_spectra = []
    mask = []
    
    for i in range(1, len(data_lattice) - 1):  # Пропускаем границы
        # Определяем границы интервала
        left_bound = w_data_midpoints[i - 1]
        right_bound = w_data_midpoints[i]
        
        # Находим точки из model_lattice, которые попадают в этот интервал
        indices_in_interval = np.where((model_lattice >= left_bound) & (model_lattice <= right_bound))[0]
        points_in_interval = model_lattice[indices_in_interval]
        intensities_in_interval = spectra[indices_in_interval]
        
        if len(points_in_interval) > 1:  # Только если больше одной точки, можно интегрировать
            # Интегрируем по правилу трапеций
            integrated_value = np.trapz(intensities_in_interval, points_in_interval)/abs(right_bound-left_bound)
            
            mask.append(i)
            w_data_result.append(data_lattice[i])
            integrated_spectra.append(integrated_value)
    
    return np.array(w_data_result), np.array(integrated_spectra), np.array(mask)

def mt_ckd_h2o_absco(p_atm, t_atm, h2o_vmr, wv1abs, wv2abs, dvabs, dat, frgnx='1', radflag=True):

    def calculate_radiation_term(xvi, temp):

        xkt = temp / RADCN2
        nvi = len(xvi)
        
        # Локальные переменные
    
        # xvi = np.array(vi)
        # print(nvi)
        # print(xvi)
        rad = np.copy(xvi)
        xviokt = xvi / xkt
    
        # Вычисления на основе условий
        rad[xviokt <= 0.01] = 0.5 * xviokt[xviokt <= 0.01] * xvi[xviokt <= 0.01]
        mask_elsewhere = (xviokt > 0.01) & (xviokt <= 10.0)
        expvkt = np.exp(-xviokt[mask_elsewhere])
        rad[mask_elsewhere] = xvi[mask_elsewhere] * (1 - expvkt) / (1 + expvkt)
        # print(rad)
        # print(xviokt)
        return rad
    
    ref_press = dat.variables["ref_press"][:]
    ref_temp = dat.variables["ref_temp"][:]
    wavenumber = dat.variables["wavenumbers"][:]
    self_absco_ref = dat.variables["self_absco_ref"][:]
    self_texp = dat.variables["self_texp"][:]
    # print(ref_press)
    if frgnx == '1':
        for_absco_ref = dat.variables["for_closure_absco_ref"][:]
    else:
        for_absco_ref = dat.variables["for_absco_ref"][:]

    nwv = int((wv2abs - wv1abs) / dvabs) + 1

    # Шаг спектра
    dvc = wavenumber[1] - wavenumber[0]

    # Определение диапазона коэффициентовоправд
    i1 = np.searchsorted(wavenumber, wv1abs - 2 * dvc, side='left')
    i2 = np.searchsorted(wavenumber, wv2abs + 2 * dvc, side='right')
    ncoeff = i2 - i1

    # Вычисление коэффициентов
    rho_rat = (p_atm / ref_press) * (ref_temp / t_atm)
    sh2o_coeff = self_absco_ref[i1:i2] * (ref_temp / t_atm) ** self_texp[i1:i2]
    sh2o_coeff *= h2o_vmr * rho_rat

    if radflag:
        rad = calculate_radiation_term(wavenumber[i1:i2], t_atm)
        sh2o_coeff *= rad

    # Интерполяция на выходную сетку
    nptabs = int((wv2abs - wv1abs) / dvabs + 1)
    self_absco = np.interp(
        np.linspace(wv1abs, wv2abs, nptabs),
        wavenumber[i1:i2],
        sh2o_coeff,
    )

    fh2o_coeff = for_absco_ref[i1:i2] * (1 - h2o_vmr) * rho_rat
    if radflag:
        fh2o_coeff *= rad

    for_absco = np.interp(
        np.linspace(wv1abs, wv2abs, nptabs),
        wavenumber[i1:i2],
        fh2o_coeff,
    )

    return self_absco, for_absco



