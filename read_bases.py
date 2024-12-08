import os
import numpy as np
import zipfile

## ======================================================
def read_hitran2012(filename):

    if not os.path.exists:
        raise ImportError('The input filename"' + filename + '" does not exist.')

    if filename.endswith('.zip'):
        zip = zipfile.ZipFile(filename, 'r')
        (object_name, ext) = os.path.splitext(os.path.basename(filename))
        filehandle = zip.read(object_name).splitlines()
    else:
        filehandle = open(filename, 'r')

    data = {'M':[],               ## molecule identification number
            'I':[],               ## isotope number
            'linecenter':[],      ## line center wavenumber (in cm^{-1})
            'S':[],               ## line strength, in cm^{-1} / (molecule m^{-2})
            'Acoeff':[],          ## Einstein A coefficient (in s^{-1})
            'gamma-air':[],       ## line HWHM for air-broadening
            'gamma-self':[],      ## line HWHM for self-emission-broadening
            'Epp':[],             ## energy of lower transition level (in cm^{-1})
            'N':[],               ## temperature-dependent exponent for "gamma-air"
            'delta':[],           ## air-pressure shift, in cm^{-1} / atm
            'Vp':[],              ## upper-state "global" quanta index
            'Vpp':[],             ## lower-state "global" quanta index
            'Qp':[],              ## upper-state "local" quanta index
            'Qpp':[],             ## lower-state "local" quanta index
            'Ierr':[],            ## uncertainty indices
            'Iref':[],            ## reference indices
            'flag':[],            ## flag
            'gp':[],              ## statistical weight of the upper state
            'gpp':[]}             ## statistical weight of the lower state

    print('Reading "' + filename + '" ...')
    
## ======================================================
    for line in filehandle:
        if (len(line) < 160):
            raise ImportError('The imported file ("' + filename + '") does not appear to be a HITRAN2012-format data file.')

        data['M'].append(int(line[0:2]))
        data['I'].append(int(line[2]))
        data['linecenter'].append(float(line[3:15]))
        data['S'].append(float(line[15:25]))
        data['Acoeff'].append(float(line[25:35]))
        data['gamma-air'].append(float(line[35:40]))
        data['gamma-self'].append(float(line[40:45]))
        data['Epp'].append(float(line[45:55]))
        data['N'].append(float(line[55:59]))
        data['delta'].append(float(line[59:67]))
        data['Vp'].append(line[67:82])
        data['Vpp'].append(line[82:97])
        data['Qp'].append(line[97:112])
        data['Qpp'].append(line[112:127])
        data['Ierr'].append(line[127:133])
        data['Iref'].append(line[133:145])
        data['flag'].append(line[145])
        data['gp'].append(line[146:153])
        data['gpp'].append(line[153:160])

    if filename.endswith('.zip'):
        zip.close()
    else:
        filehandle.close()

    for key in data:
        data[key] = np.array(data[key])

    return(data)





def read_bezard(filename):

    if not os.path.exists:
        raise ImportError('The input filename"' + filename + '" does not exist.')

    if filename.endswith('.zip'):
        zip = zipfile.ZipFile(filename, 'r')
        (object_name, ext) = os.path.splitext(os.path.basename(filename))
        filehandle = zip.read(object_name).splitlines()
    else:
        filehandle = open(filename, 'r')

    data = {'I':[],               ## isotope number
            'linecenter':[],      ## line center wavenumber (in cm^{-1})
            'S':[],               ## line strength, in cm^{-1} / (molecule m^{-2})
            'gamma-air':[],       ## line HWHM for air-broadening
            'gamma-self':[],      ## line HWHM for self-emission-broadening
            'Epp':[],             ## energy of lower transition level (in cm^{-1})
            'N':[],               ## temperature-dependent exponent for "gamma-air"
            'delta':[]}           ## air-pressure shift, in cm^{-1} / atm

    print('Reading "' + filename + '" ...')

    for line in filehandle:

        data['linecenter'].append(float(line[1:10]))
        data['S'].append(float(line[11:20]))
        data['Epp'].append(float(line[26:35]))
        data['I'].append(int(line[76:79]))

    if filename.endswith('.zip'):
        zip.close()
    else:
        filehandle.close()

    for key in data:
        data[key] = np.array(data[key])

    return(data)




## ======================================================
def read_ames(filename):

    if not os.path.exists:
        raise ImportError('The input filename"' + filename + '" does not exist.')

    if filename.endswith('.zip'):
        zip = zipfile.ZipFile(filename, 'r')
        (object_name, ext) = os.path.splitext(os.path.basename(filename))
        filehandle = zip.read(object_name).splitlines()
    else:
        filehandle = open(filename, 'r')

    data = {'I':[],               ## isotope number
            'linecenter':[],      ## line center wavenumber (in cm^{-1})
            'S':[],               ## line strength, in cm^{-1} / (molecule m^{-2})
            'Acoeff':[],          ## Einstein A coefficient (in s^{-1})
            'gamma-self':[],      ## line HWHM for self-emission-broadening
            'gamma-air':[],       ## line HWHM for air-broadening
            'Epp':[],             ## energy of lower transition level (in cm^{-1})
            'N':[],               ## temperature-dependent exponent for "gamma-air"
            'delta':[]}           ## air-pressure shift, in cm^{-1} / atm

    print('Reading "' + filename + '" ...')

    for line in filehandle:

        data['linecenter'].append(float(line[2:15]))
        data['S'].append(float(line[15:26]))
        data['Acoeff'].append(float(line[26:37]))
        data['Epp'].append(float(line[54:66]))
        data['I'].append(int(line[0:2]))
        data['gamma-self'].append(float(line[123:128]))
        data['gamma-air'].append(float(line[37:42]))
        data['N'].append(float(line[128:131]))
        # data['N'].append(float(line[42:46]))
        data['delta'].append(float(line[46:54]))

    if filename.endswith('.zip'):
        zip.close()
    else:
        filehandle.close()

    for key in data:
        data[key] = np.array(data[key])

    return(data)

## ======================================================

def read_hot_ames(filename):

    if not os.path.exists(filename):
        raise ImportError('The input filename "' + filename + '" does not exist.')

    if filename.endswith('.zip'):
        zip = zipfile.ZipFile(filename, 'r')
        object_name, ext = os.path.splitext(os.path.basename(filename))
        filehandle = zip.read(object_name).splitlines()
    else:
        filehandle = open(filename, 'r')

    data = {
        'I': [],            # isotope number
        'linecenter': [],   # line center wavenumber (in cm^{-1})
        'S': [],            # line strength
        'gamma-self': [],   # line HWHM for self-emission-broadening
        'gamma-air': [],       ## line HWHM for air-broadening
        'N': [],            # temperature-dependent exponent for self-broadening
        'Epp': [],          # energy of lower transition level (in cm^{-1})
        'delta': []         # air-pressure shift (delta)
    }

    print('Reading "' + filename + '" ...')

    for line in filehandle:
        if not line:
            continue
        
        data['I'].append(int(line[1:3]))
        data['linecenter'].append(float(line[3:16]))
        data['S'].append(float(line[16:27]))
        data['gamma-self'].append(float(line[125:130]))  # gam self at 350 K
        data['gamma-air'].append(float(line[37:42]))
        data['N'].append(float(line[130:133]))    # n self at 350 K
        data['Epp'].append(float(line[55:67]))
        data['delta'].append(float(line[47:55]))

    if filename.endswith('.zip'):
        zip.close()
    else:
        filehandle.close()

    for key in data:
        data[key] = np.array(data[key])

    return data
