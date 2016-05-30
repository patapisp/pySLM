__author__ = 'Chronis'
import numpy as np
from scipy.special import j1
from astropy.io import fits
import pygame, os, time, pickle
import win32api
from tkinter import *
from tkinter import ttk
from tkinter import filedialog
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.image as mplimg
from matplotlib.colors import Normalize
from matplotlib import cm
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
matplotlib.use('TkAgg')


class SLM:
    """
    A class that enables to write an image to the Holoeye SLM seen as a second screen
    with dimensions 1024x768
    """
    def __init__(self):
        # get path that file is saved in
        self.path = os.path.dirname(os.path.realpath("__file__"))
        # find connected screens
        screens = win32api.EnumDisplayMonitors()
        if len(screens) < 2:
            raise UserWarning('No second screen connected')
        self.hdc = screens[1][0]  # handle of second screen PYHANDLE
        self.dim = screens[1][2]  # screen dimensions
        self.left = self.dim[0]
        self.top = self.dim[1]
        self.right = self.dim[2]
        self.bottom = self.dim[3]
        self.width = abs(self.right-self.left)
        self.height = abs(self.bottom-self.top)
        self.size = (self.width, self.height)
        self.dimensions = (self.width, self.height, 3)
        # set Windows ENVIRONMENTAL VARIABLE of SDL window position to the top-left corner
        # of the second screen
        os.environ['SDL_VIDEO_WINDOW_POS'] = "%d,%d" % (self.left, self.top)
        pygame.display.init()
        self.screen = pygame.display.set_mode(self.size)
        # create surface object
        self.surface = pygame.surfarray.make_surface(np.zeros(self.dimensions, dtype=np.uint8))
        self.screen.blit(self.surface, (0, 0))
        pygame.display.flip()
        self.pxarray = []
        self.maps = self.load_maps()  # dictionary to store phase maps imported from file
        self.active = False

        # SLM pixel size
        self.pixelSize = 36

    def create_reference_array(self):
        """
        Reference surface pixels with numpy array. Every change made in the array will
        automatically change the pixels in the surface. Ideal for fast operations.
        :return:
        """
        self.active = True
        self.pxarray = pygame.surfarray.pixels3d(self.surface)
        pygame.display.flip()
        return

    def delete_reference_array(self):
        """
        Delete previously referenced array
        :return:
        """
        del self.pxarray
        return

    def load_maps(self):
        """
        Import previously used phase maps that are stored in phase_maps.p
        :return:
        """
        fname = self.path + '\\phase_maps.p'
        if os.path.isfile(fname):
            return pickle.load(open(fname, 'rb'))
        else:
            return {}

    def save_maps(self):
        """
        Save imported map to phase_maps dictionary pickle file
        :return:
        """
        fname = self.path + '\\phase_maps.p'
        pickle.dump(self.maps, open(fname, 'wb'))
        return

    def update(self):
        """
        Get buffer of second screen
        :return:
        """
        pygame.display.update()
        return

    def quit(self):
        """
        Quits display
        :return:
        """
        pygame.display.quit()
        return

    def draw(self, p):
        """
        Draw array onto second screen.
        :param p: p is an uint8 numpy array with the correct dimensions
        :return:
        """
        self.active = True
        surface = pygame.surfarray.make_surface(p)
        self.screen.blit(surface, (0, 0))
        pygame.display.flip()
        return

    def four_quadrants(self, value, center):
        """
        Make a matrix representing the four quadrants phase map
        with four 0,value alternating regions
        :param value: Usually value[0]=255, value[1]=pi but can change
        :param center: Center (x0,y0) of four quadrants
        :return:None
        """
        if not((0 < value[0] < 256) and (0 < value[1] < 256)):
            raise ValueError('value must be unsigned 8bit not %i' % value)
        _x = int(center[0])
        _y = int(center[1])
        four_q = np.zeros((self.width, self.height), dtype=np.uint8)*value[0]
        four_q[_x:-1, 0:_y] = value[1]
        four_q[0:_x, _y:-1] = value[1]
        return four_q

    def import_phase_map(self, file):
        """
        Imports an phase map stored in a file and saves it in an
        numpy array for use on the SLM
        :param file: A .txt or .fits file of a phase map
        :return: write an 8bit array of the correct dimensions in maps dict
        """
        filename, extension = os.path.splitext(file)
        if extension == '.txt':
            p = np.loadtxt(file, dtype=int, delimiter=' ')
            if np.shape(p) != (1024, 768):
                return
            if np.shape(p) == (768, 1024):
                p = p.T
            m = np.zeros((1024, 768, 3), dtype=np.uint8)
            m[:, :, 0] = p
            m[:, :, 1] = p
            m[:, :, 2] = p
            self.maps[os.path.basename(filename)] = {'data': m}
        elif extension == '.fits':
            hdu = fits.open(file)
            p = hdu[0].data
            m = np.zeros((1024, 768, 3), dtype=np.uint8)
            m[:, :, 0] = p
            m[:, :, 1] = p
            m[:, :, 2] = p
            self.maps[os.path.basename(filename)] = {'data': m}
        elif extension == '.bmp':
            p = pygame.image.load(file)
            p = pygame.surfarray.array3d(p)
            self.maps[os.path.basename(filename)] = {'data': p}
        else:
            raise UserWarning('File extension %s not supported.' % extension)

        # if np.shape(self.maps[filename]) != (1024, 768, 3):
        #    raise UserWarning('Wrong dimensions for SLM of size %i x %i (x %i)' % (self.width, self.height, 3))

        self.save_maps()
        return filename

    def airy(self, r, F, l, I0):
        """
        Airy pattern
        :param r: radial distance from star
        :param F: F number
        :param l: wavelength
        :param I0: Max intensity
        :return:
        """
        if r == 0:
            return I0
        x = np.pi*r*self.pixelSize/(l*1e-3*F)
        I = I0*(2*j1(x)/x)**2
        return I

    def lorentzian(self, r, a, b):
        """
        Creates the value of the lorentzian at distance r with center at a and scale b
        and re-normalizes to 1
        :param r: distance from star center
        :param a: position of star center
        :param b: scale/intensity
        :return: value
        """
        l = 1/(np.pi*b*(((r-a)**2/b**2)+1))
        norm = 1/(np.pi*b)
        l /= norm
        return l

    def four_qs(self, xp, yp, c, val1, val2):
        """
        Calculates abstract four quadrants value for pixel x,y
        :param xp: x in pixel coordinates
        :param yp: y in pixel coordinates
        :param c: (x,y) of star center
        :param val1: value 1 of 4 quadrants phase mask
        :param val2: value 2 of 4 quadrants phase mask
        :return:
        """
        xs, ys = c
        x = xs-xp
        y = ys-yp
        if (x == 0 and y >= 0) or (y==0 and x > 0):
            return val2
        elif x == 0 and y<0 or (y==0 and x <= 0):
            return val1
        else:
            pass
        phi = np.arctan(y/x)
        expr = (0 <= phi <= 0.5*np.pi) or (np.pi < phi < 1.5*np.pi)
        if expr:  # pixel is in first or third quadrant
            return val2
        else:  # pixel is in second or fourth quadrant
            return val1
        
    def eight_octants(self, xp, yp, c, val1, val2):
        """
        Method that creates 8octants phase mask
        :param xp: pixel x coordinate
        :param yp: pixel y coordinate
        :param c: center coordinates in tuple (xc, yc)
        :param val1: gray value 1
        :param val2: gray value 2
        :return:
        """
        xs, ys = c
        x = xp-xs
        y = yp-ys

        expr2 = lambda phi:(0.25*np.pi < phi)
        expr1 = lambda phi:(0.25*np.pi <= phi)
        expr3 = lambda phi:(-0.5*np.pi < phi < -0.25*np.pi)
        expr4 = lambda phi:(-0.25*np.pi >= phi > -0.5*np.pi)
        phi = 0

        if x == 0 and y < 0:
            return val2
        elif x == 0 and y >= 0:
            return val1
        elif y==0 and x<0:
            return val2

        else:
            phi = np.arctan(y/x)
            if y > 0 and x > 0:
                if expr2(phi):  # pixel is in first or third quadrant
                    return val1
                else:  # pixel is in second or fourth quadrant
                    return val2
            elif y < 0 and x < 0:
                if expr1(phi):  # pixel is in first or third quadrant
                    return val1
                else:  # pixel is in second or fourth quadrant
                    return val2
            elif y < 0 and x > 0:
                if expr3(phi):  # pixel is in first or third quadrant
                    return val2
                else:  # pixel is in second or fourth quadrant
                    return val1
            else:
                if expr4(phi):  # pixel is in first or third quadrant
                    return val2
                else:  # pixel is in second or fourth quadrant
                    return val1

    def eight_octants_old(self, xp, yp, c, val1, val2):
        """

        :param xp: pixel x coordinate
        :param yp: pixel y coordinate
        :param c: center coordinates in tuple (xc, yc)
        :param val1: gray value 1
        :param val2: gray value 2
        :return:
        """
        xs, ys = c
        x = xp-xs
        y = yp-ys

        expr2 = lambda phi:(0.25*np.pi < phi)
        expr1 = lambda phi:(0.25*np.pi <= phi)
        expr3 = lambda phi:(-0.5*np.pi < phi < -0.25*np.pi)
        expr4 = lambda phi:(-0.75 > phi > -1*np.pi)
        phi = 0
        if (x == -1 and y == 0) or (x == 0 and y == -1):
            return val1
        elif (x == 0 and y == 0) or (x == -1 and y == -1):
            return val2
        if x == 0 and y < 0:
            return val2
        elif x == 0 and y >= 0:
            return val1

        else:
            phi = np.arctan(y/x)
            if y > 0 and x > 0:
                if expr2(phi):  # pixel is in first or third quadrant
                    return val1
                else:  # pixel is in second or fourth quadrant
                    return val2
            elif y < 0 and x < 0:
                if expr1(phi):  # pixel is in first or third quadrant
                    return val1
                else:  # pixel is in second or fourth quadrant
                    return val2
            elif y < 0 and x > 0:
                if expr3(phi):  # pixel is in first or third quadrant
                    return val2
                else:  # pixel is in second or fourth quadrant
                    return val1
            else:
                if expr4(phi):  # pixel is in first or third quadrant
                    return val2
                else:  # pixel is in second or fourth quadrant
                    return val1

    def pixel_value(self, x, y, c1, c2, i1, i2, val1, val2, F1, F2, l1, l2, mask='FQPM'):
        """
        Calculates value of pixel for 4 quadrants function
        :param x: coordinates of pixel
        :param y:
        :param type: what type of map to apply
        :return:
        """
        x1, y1 = c1
        a1 = np.sqrt(x1**2+y1**2)
        x2, y2 = c2
        a2 = np.sqrt(x2**2+y2**2)
        r1 = np.sqrt((x1-x)**2 + (y1-y)**2)  # doesn't have to be an integer
        r2 = np.sqrt((x2-x)**2 + (y2-y)**2)
        k1_airy = self.airy(r1, F1, l1, i1)
        k2_airy = self.airy(r2, F2, l2, i2)
        norm_airy = k1_airy + k2_airy
        k1_airy /= norm_airy
        k2_airy /= norm_airy
        if mask == 'FQPM':
            val_airy = k1_airy*self.four_qs(x, y, c1, val1, val2) + \
                       k2_airy*self.four_qs(x, y, c2, val1, val2)
            return val_airy
        elif mask == 'EOPM':
            val_airy = k1_airy*self.eight_octants(x, y, c1, val1, val2) + \
                       k2_airy*self.eight_octants(x, y, c2, val2, val1)
            return val_airy

    def open_real_psf(self, psf_file):
        # choose image of real psf
        hdu = fits.open(psf_file)
        self.psf = hdu[0].data
        self.psf = self.psf.T
        try:
            c = hdu[0].header['center']  # center of weighting funnction
            c = c.split(sep=',')
            y, x = c
            self.real_center = [int(x), int(y)]
        except:
            print("center not found in header")
            return None
        return

    def real_psf_weight(self, x, y, xc, yc):
        """
        Find out how much pixel is from center and get value of weight from real psf
        :param x:
        :param y:
        :param c: center of star to get weight for
        :return: weight
        """
        x_i = x - xc + self.real_center[0]
        y_i = y - yc + self.real_center[1]
        return self.psf[x_i, y_i]



    def pixel_value_real_psf(self, x, y, c1, c2, val1, val2, i1, i2, psf_file):
        """
        Value of pixel in phase mask for binary case
        :param x:
        :param y:
        :param c1: center of primary star
        :param c2: center of second star
        :param val1: phase value 1
        :param val2: phase value 2
        :param i1: intensity of primary star
        :param i2: intensity of secondary star
        :param psf_file: file containing the already normalized and smoothed real psf
        :return: value of phase mask
        """
        x1, y1 = c1
        x2, y2 = c2
        r1 = np.sqrt((x1-x)**2 + (y1-y)**2)  # doesn't have to be an integer
        r2 = np.sqrt((x2-x)**2 + (y2-y)**2)
        k1_real = i1*self.real_psf_weight(x, y, x1, y1)
        k2_real = i2*self.real_psf_weight(x, y, x2, y2)
        norm_real = k1_real + k2_real
        k1_real /= norm_real
        k2_real /= norm_real
        val_airy = k1_real*self.four_qs(x, y, c1, val1, val2) + k2_real*self.four_qs(x, y, c2, val1, val2)
        return val_airy



















