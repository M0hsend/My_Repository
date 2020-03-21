import numpy as np
from scipy import constants as pc
from numpy.fft import fftshift, ifft2, ifftshift
import matplotlib.pyplot as plt
import h5py


def rad_array(arr_len):
    """
    create array with values equal to distance in px from centre
    :param arr_len:
    :return:
    """
    half_len = arr_len / 2
    v = np.linspace(-half_len, half_len, arr_len)
    v = np.repeat(v, arr_len, axis = 0)
    v = np.reshape(v, (arr_len, arr_len))
    v = np.sqrt(v**2 + v.T**2)
    return v


def e_lambda(e_0):
    """
    relativistic electron wavelength

    :param e_0: int
        accelerating voltage in volts
    :return:
    e_lambda: float
        wavelength in meters
    """
    import numpy as np
    
    # m = 9.109383*10**-31 # electron rest mass in kg
    # e = 1.602177*10**-19 # primary charge in C
    # h = 6.62607*10**-34 # planck const in m**2*kg/s
    # c =  299792458 # speed of light in m/s
    
    e_lambda = (pc.h * pc.c) / np.sqrt((pc.e * V)**2  + 2 * pc.e * V * pc.m_e * pc.c**2)
    
    return e_lambda


def FuncAberrUV(u, v, aberrcoeff):
    """
    Returns the aberration function using the aberration coefficients
    :param u: fourier coordinate x
    :param v: fourier coordinate y
    :param aberrcoeff: array of floats
        array of aberration coefficients
    :return:
    fun_aberr:
    """

    u2 = u*u
    u3 = u2*u
    u4 = u3*u
    
    v2 = v*v
    v3 = v2*v
    v4 = v3*v
    
    # aberr are in unit of meter.
    C1   = aberrcoeff[0] # defocus
    C12a = aberrcoeff[1] # 2 stig
    C12b = aberrcoeff[2] # 2 stig
    C23a = aberrcoeff[3] # 3 stig
    C23b = aberrcoeff[4] # 3 stig
    C21a = aberrcoeff[5] # coma 
    C21b = aberrcoeff[6] # coma
    C3   = aberrcoeff[7] # Spherical abb
    C34a = aberrcoeff[8] # 4 stig
    C34b = aberrcoeff[9] # 4 stig
    C32a = aberrcoeff[10] # star
    C32b = aberrcoeff[11] # star
    
    # output:  chi function. in unit of meter*radian.  multiply by 2pi/lambda to get dimensionless
    func_aberr =  1/2*C1*(u2+v2)\
            + 1/2*(C12a*(u2-v2) + 2*C12b*u*v)\
            + 1/3*(C23a*(u3-3*u*v2) + C23b*(3*u2*v - v3))\
            + 1/3*(C21a*(u3+u*v2) + C21b*(v3+u2*v))\
            + 1/4* C3*(u4+v4+2*u2*v2)\
            + 1/4* C34a*(u4-6*u2*v2+v4)\
            + 1/4* C34b*(4*u3*v-4*u*v3)\
            + 1/4* C32a*(u4-v4)\
            + 1/4* C32b*(2*u3*v + 2*u*v3)\
    
    return func_aberr


def define_probe_function(V, alpha, px_cal, array_px, aberr_input, dose = 1, plot_me = False):
    """

    :param V: int
        accelerating voltage - volts
    :param alpha:
    :param px_cal: float
        pixel size in m
    :param array_px: int
        number of pixels in detector plane
    :param aberr_input: numpy array of floats
        aberration coeffs
    :param dose: 
    :param plot_me: boolean
        default False
    :return:
    func_probe: complex numpy array
        probe function
    """
    l = e_lambda(V)
    #% probe function
    #% chi function
    cen = array_px / 2 # center of array
    #K_max = ((cen * px_cal) - (px_cal/2)) / l # max scattering vector 
    K_px =  l / (px_cal * array_px)
    K_max = K_px * array_px
    #print(K_px, K_max)
    Kx = np.linspace(-K_max , K_max, array_px)
    Kx = np.repeat(Kx,array_px,  axis = 0)
    Kx = np.reshape(Kx, (array_px, array_px))
    Ky = np.copy(Kx)
    Kx = Kx.T
    func_aber = FuncAberrUV(Kx,Ky,aberr_input)
    
    # transfer function
    func_transfer=np.exp((-1j*2*np.pi/ (l)) * func_aber)
    
    # aperture function
    #func_ObjApt = ones(size(ptycho.Kwp));
    #func_ObjApt( ptycho.Kwp > ptycho.ObjApt_angle) = 0;
    #array_px = Kx.shape
    func_ObjApt = np.zeros((array_px, array_px), dtype = int)
    xx, yy = np.mgrid[:array_px, :array_px]
    
    circle =np.sqrt((xx - cen) ** 2 + (yy - cen) ** 2)
    alpha_px = alpha / K_px # alpha in px
    func_ObjApt[np.where(circle< alpha_px)] = 1
    #% dose equals to the summed intensity of the average ronchigram.
    #dose = sum(ptycho.pacbed(:)) *1.0;
    #% for resampled ronchigram
    #% pacbed_rs = interp2(tx_wp,ty_wp,mean_m_wp,Kx,Ky,'cubic',0);
    #% dose = sum(pacbed_rs(:)) *1.0;
    
    #% normalize the Objective Aperture to the desired number of electrons
    scaling = np.sqrt(dose/func_ObjApt.sum())
    func_ObjApt = func_ObjApt* scaling
    
    #% convergence aperture function; to filter out the updated probe func
    #% during ptychography iterations. 
    #% ObjAptMask = func_ObjApt./sum(func_ObjApt(:));
    
    #% probe function - reciprocal space
    A = func_ObjApt*func_transfer
    #% probe function - real space
    func_probe = fftshift(ifft2(ifftshift(A)))
    
    if plot_me:
        fig_mul = 0.0625
        #max_fig = px_cal * array_px
        im_lim = [cen - (array_px * fig_mul), cen + (array_px * fig_mul)]
        fig_lim =[-px_cal *array_px *fig_mul ,px_cal *array_px *fig_mul]
        plt.figure()
        fig, [[ax11, ax12], [ax21, ax22]] = plt.subplots(2, 2)
        ax11.imshow(np.angle(A))
        ax11.set_title('Aperture Phase Surface')
        ax11.set_xlim(im_lim)
        ax11.set_ylim(im_lim)
        ax12.imshow(abs(func_probe))
        ax12.set_title('Probe Function')
        ax12.set_xlim(im_lim)
        ax12.set_ylim(im_lim)
        ax21.imshow(np.real(func_probe))
        ax21.set_title('real')
        ax21.set_xlim(im_lim)
        ax21.set_ylim(im_lim)
        ax22.imshow(np.imag(func_probe))
        ax22.set_title('imag')
        ax22.set_xlim(im_lim)
        ax22.set_ylim(im_lim)

        x_list = np.arange(-cen * px_cal, cen*px_cal, px_cal)
        plt.figure()
        plt.plot(x_list, abs(func_probe[int(array_px/2), :]))
        plt.plot(x_list, np.real(func_probe[int(array_px/2), :]))
        plt.plot(x_list, np.imag(func_probe[int(array_px/2), :]))
        plt.xlim(fig_lim)
        plt.title('df = ' + str(aberr_input[0]) + ' Cs = ' + str(aberr_input[7]))
    return func_probe   


#%%
V = 80000
alpha = 10e-3#100e-3 # rad
px_cal = 0.2e-10#0.45e-10 # in m 
array_px = 72
output_size = 256
save_hdf5 = False
save_path = r'Y:\2019\cm22979-8\processing\Merlin\20191114_15kVptycho_graphene\probe_sims'
save_file = r'\15kV_10um_Cs987um'

aberrcoeff = np.zeros((12))
aberrcoeff[0] = 200e-10 # defocus
aberrcoeff[1] = 0# 2 stig
aberrcoeff[2] = 0  # 2 stig
aberrcoeff[3] = 0 # 3 stig
aberrcoeff[4] = 0 # 3 stig
aberrcoeff[5] = 0 # coma 
aberrcoeff[6] = 0 # coma
aberrcoeff[7] = 0#987e-6 # Spherical abb
aberrcoeff[8] = 0# 4 stig
aberrcoeff[9]  = 0# 4 stig
aberrcoeff[10] = 0 # star
aberrcoeff[11]  = 0 # star
aberr_input = aberrcoeff
probe = define_probe_function(V,alpha, px_cal,array_px, aberrcoeff, dose = 1, plot_me = True)
px_from, px_to = int(array_px/2 - output_size / 2) , int(array_px/2 + output_size / 2)
output_probe = probe#[px_from:px_to, px_from:px_to]
plt.figure(); plt.imshow(np.real(output_probe))

if save_hdf5 == True:
    output_probe = output_probe[np.newaxis, np.newaxis, np.newaxis, np.newaxis,np.newaxis, :,:] 
    fn = save_path + save_file
    d5 = h5py.File(fn +'.hdf5' , 'w')
    d5.create_dataset('entry_1/process_1/output_1/probe', data = output_probe)
    d5.create_dataset('entry_1/process_1/PIE_1/detector/binning', data = [1,1])
    d5.create_dataset('entry_1/process_1/PIE_1/detector/upsample', data = [1,1])
    d5.create_dataset('entry_1/process_1/PIE_1/detector/crop',data = [output_size, output_size])
    d5.create_dataset('entry_1/process_1/common_1/dx', data = [4.52391605e-11 , 4.52391605e-11 ]) # px size
    d5.close()