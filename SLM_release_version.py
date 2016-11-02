__author__ = 'Chronis'
from pySLM.definitions import SLM
import numpy as np
from tkinter import _setit
import PIL
from astropy.io import fits
import pygame, os, time, pickle
from tkinter import *
from tkinter import ttk
from tkinter import filedialog
import matplotlib
import matplotlib.pyplot as plt
import threading
import matplotlib.image as mplimg
from matplotlib.colors import Normalize
from matplotlib import cm
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
matplotlib.use('TkAgg')


def cart2pol(x,y):
    """
    Takes cartesian (2D) coordinates and transforms them into polar.
    """
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    return (rho, phi)


class dummyClass:
    def __init__(self):
        print('Dummy class')
        self.maps = {'zero': np.zeros((1024, 768, 3))}
        self.SLM_type = 'None'
        self.pixelSize = 8
        self.dimensions = (1024, 768, 3)
        self.width = 1024
        self.height = 768
        self.size = (1024, 768)


class StdoutRedirector(object):
    def __init__(self, text_widget):
        self.text_space = text_widget

    def write(self, string):
        self.text_space.insert('end', string)
        self.text_space.see('end')

    def flush(self):
        pass


def array2PIL(arr, size):
    mode = 'RGBA'
    arr = arr.reshape(arr.shape[0]*arr.shape[1], arr.shape[2])
    if len(arr[0]) == 3:
        arr = np.c_[arr, 255*np.ones((len(arr), 1), np.uint8)]
    return PIL.Image.frombuffer(mode, size, arr.tostring(), 'raw', mode, 0, 1)


class DropMenu:
    """
    DropMenu is a widget that will contain various functionalities of a menu
    """
    def __init__(self, master, window):
        # Create dropdown menu
        self.path = os.getcwd()
        self.window = window
        self.master = master
        self.menu = Menu(self.master)
        self.master.config(menu=self.menu)
        # File Option************************************************
        self.FileMenu = Menu(self.menu)
        self.menu.add_cascade(label='File', menu=self.FileMenu)

        self.FileMenu.add_command(label='Open phase map')
        self.FileMenu.add_command(label='Save as FITS', command=lambda: self.save_fits())
        self.FileMenu.add_command(label='Save weighting function')
        self.FileMenu.add_separator()
        self.FileMenu.add_command(label='Exit', command=self._quit)
        # Settings option***********************************************
        self.SettingsMenu = Menu(self.menu)
        self.menu.add_cascade(label='Settings', menu=self.SettingsMenu)

        self.SettingsMenu.add_command(label='Calibration curve', command=self.calibration_callback)
        self.SettingsMenu.add_command(label='Star info', command=self.star_info_callback)
        # Tools option**************************************************
        self.ToolMenu = Menu(self.menu)
        self.menu.add_cascade(label='Tools', menu=self.ToolMenu)

        self.ToolMenu.add_command(label='Count')
        self.ToolMenu.add_command(label='Histogram')
        # Help option ********************************************
        self.HelpMenu = Menu(self.menu)
        self.menu.add_cascade(label='Help', menu=self.HelpMenu)

        self.HelpMenu.add_command(label='Documentation')
        self.HelpMenu.add_command(label='App Help')

        # Variables **********************************************
        try:
            self.menu_data = pickle.load(open("SLM_data.p", 'rb'))
            self.phase_curve = self.menu_data['phase curve']
        except:
            file = filedialog.askopenfilename(title="Select phase curve(.npy)")
            phase = np.load(file)
            self.menu_data = {'phase curve': phase}
            self.phase_curve = phase
            pickle.dump(self.menu_data, open("SLM_data.p", 'wb'))

        p = np.polyfit(self.phase_curve, np.arange(0, 256), deg=5)
        self.rad_2_gray = np.poly1d(p)
        self.v0, self.vpi2, self.vpi, self.v3pi2 = StringVar(), StringVar(), StringVar(), StringVar()
        self.v0.set('0')
        self.vpi2.set('111')
        self.vpi.set('200')
        self.v3pi2.set('255')
        self.slm_pxl = StringVar()
        self.slm_pxl.set('36')
        self.intensity = StringVar()
        self.wavelength = StringVar()
        self.Fnum = StringVar()
        self.lD = StringVar()
        self.lD.set('4')

    def star_info_callback(self):
        toplevel_r = Toplevel()
        toplevel_r.title('Star info')
        toplevel_r.geometry("400x150+300+300")
        toplevel = ttk.Frame(toplevel_r)
        toplevel.grid(column=0, row=0, sticky=(N, W, E, S))
        self.wavelength.set('633')
        wavelength_entry = Entry(toplevel, textvariable=self.wavelength,justify='center')
        wavelength_lab = Label(toplevel, text='Wavelength (nm):')
        self.Fnum.set('230')
        Fnum_entry = Entry(toplevel, textvariable=self.Fnum, justify='center')
        Fnum_lab = Label(toplevel, text='F # :')
        self.intensity.set('1')
        intensity_entry = Entry(toplevel, textvariable=self.intensity, justify='center')
        intensity_lab = Label(toplevel, text='Intensity :')
        """As discussed, here are the correct parameters for the coordinates conversion in the SLM plane :
        F# = 230
        Pixel_size = 36 um
        The spot size in the SLM plane right now is lambda*F# ~ 145 um ~ 4 pixels.
        """
        slm_pxl_lab = Label(toplevel, text='SLM pixel size (um):', justify='center')
        slm_pxl_entry = Entry(toplevel, textvariable=self.slm_pxl)
        lD_lab = Label(toplevel, text='#pixels per l/D:')
        lD_entry = Entry(toplevel, textvariable=self.lD)
        separator = ttk.Separator(toplevel, orient=VERTICAL)
        set_button = ttk.Button(toplevel, text='Set', command=self.apply_star_info)
        wavelength_lab.grid(column=0, row=0)
        wavelength_entry.grid(column=1, row=0)
        Fnum_lab.grid(column=0, row=1)
        Fnum_entry.grid(column=1, row=1)
        intensity_lab.grid(column=0, row=2)
        intensity_entry.grid(column=1, row=2)
        separator.grid(column=2, row=0, rowspan=3, sticky=(N, S))
        slm_pxl_lab.grid(column=3, row=0)
        slm_pxl_entry.grid(column=3, row=1)
        lD_lab.grid(column=3, row=2)
        lD_entry.grid(column=3, row=3)

        set_button.grid(column=0, row=4)



    def apply_star_info(self):
        pass

    def calibration_callback(self):
        toplevel_r = Toplevel()
        toplevel_r.title('Grayvalues calibration')
        toplevel_r.geometry("300x150+300+300")
        toplevel = ttk.Frame(toplevel_r)
        toplevel.grid(column=0, row=0, sticky=(N, W, E, S))

        # v0_lab = Label(toplevel, text='0')
        # v0_edit = Entry(toplevel, textvariable=self.v0, justify='center')
        # vpi2_lab = Label(toplevel, text='pi/2')
        # vpi2_edit = Entry(toplevel, textvariable=self.vpi2, justify='center')
        # vpi_lab = Label(toplevel, text='pi')
        # vpi_edit = Entry(toplevel, textvariable=self.vpi, justify='center')
        # v3pi2_lab = Label(toplevel, text='3pi/2')
        # v3pi2_edit = Entry(toplevel, textvariable=self.v3pi2, justify='center')
        # real_value_lab = Label(toplevel, text='Real values')
        # grayval_lab = Label(toplevel, text='Gray values:')
        # self.data_real = [0, 0.5*np.pi, np.pi, 1.5*np.pi]
        # self.data_gray = [0, 111, 200, 255]

        self.curve_plot, self.ax = plt.subplots(figsize=(3,3))
        self.ax.plot(np.arange(256), self.phase_curve, 'o')
        self.ax.set_xlim([-1, 260])
        self.ax.set_xlabel("gray values")
        self.ax.set_ylabel("Phase shift[$\pi$]")
        data_plot = FigureCanvasTkAgg(self.curve_plot, master=toplevel)
        data_plot.show()
        import_curve_button = ttk.Button(toplevel, text='Import curve', command=self.import_curve_callback)
        # real_value_lab.grid(column=0, row=0)
        # grayval_lab.grid(column=0, row=1)
        # v0_lab.grid(column=1, row=0)
        # v0_edit.grid(column=1, row=1)
        # vpi2_lab.grid(column=2, row=0)
        # vpi2_edit.grid(column=2, row=1)
        # vpi_lab.grid(column=3, row=0)
        # vpi_edit.grid(column=3, row=1)
        # v3pi2_lab.grid(column=4, row=0)
        # v3pi2_edit.grid(column=4, row=1)
        import_curve_button.grid(column=0, row=2)
        data_plot.get_tk_widget().grid(column=1, row=2, columnspan=4, rowspan=4)
        return

    def import_curve_callback(self):
        file = filedialog.askopenfilename(title="Select phase curve(.npy)")
        phase = np.load(file)
        self.menu_data = {'phase curve': phase}
        self.phase_curve = phase
        self.ax.set_ydata(phase)
        pickle.dump(self.menu_data, open("SLM_data.p", 'wb'))
        return


    def save_fits(self, name=None):
        """
        Save current open phase mask as a FITS file with the center information
        contained in the header
        """
        file = filedialog.asksaveasfilename(master=self.master, title='Save as..', initialdir=self.path)
        if file is None:
            return
        self.path = os.path.dirname(file)
        file += '.fits'
        # current = 0
        if name is None:
            current = self.window.maps_var.get()
        else:
            current = name

        if current == '':
            return

        mask = self.window.SLM.maps[current]['data']
        if self.window.active:
            mask = self.window.image
        hdu = fits.PrimaryHDU()
        hdu.data = mask[:, :, 0]
        hdu.header['center'] = str(self.window.center_position)
        if len(self.window.center_position) > 1:
            hdu.header['stars'] = str(self.window.multiple_star_position) + "\ ([l/D, azimuth]"

        """
        if mask['star info']:
            for k, val in mask['star info']:
                hdu.header[k] = val
        """
        hdu.header['DATE'] = time.strftime("%d/%m/%Y")
        hdu.writeto(file)
        return

    def _quit(self):
        self.window.SLM.quit()
        self.master.quit()     # stops mainloop
        self.master.destroy()
        return


class SLMViewer:
    """
    Basic GUI that enables communication with SLM , on/off switch and
    import/manipulation of phase maps
    """
    def __init__(self, root):
        self.master = Frame(root)
        self.master.grid(column=0, row=0, sticky=(N, W, E, S))
        root.title('SLM Controller')
        try:
            self.SLM = SLM()
            print("SLM type is %s"%self.SLM.SLM_type)
        except UserWarning:
            self.SLM = dummyClass()
            #raise UserWarning('No SLM connected.')

        self.menu = DropMenu(root, self)  # add drop-down menu
        #self.SLM.pixelSize = int(self.menu.slm_pxl.get())

        # =====================================================================================
        # make canvas
        self.off_image = np.zeros(self.SLM.dimensions)
        self.image = self.off_image
        self.fig, self.ax = plt.subplots()
        self.norm = Normalize(vmin=0, vmax=255)
        self.cmap = cm.gray
        self.im = plt.imshow(self.image[:, :, 0].T, cmap=self.cmap, norm=self.norm)
        self.ax.get_xaxis().set_visible(False)
        self.ax.get_yaxis().set_visible(False)
        # get image plot onto canvas and app
        self.data_plot = FigureCanvasTkAgg(self.fig, master=self.master)
        self.data_plot.get_tk_widget().configure(borderwidth=0)
        self.fig.suptitle('SLM type : %s'%self.SLM.SLM_type, fontsize=12, fontweight='bold')

        self.data_plot.show()
        self.fig.canvas.mpl_connect('button_press_event', self.click_callback)
        # ====================================================================================
        # import phase maps frame
        self.import_maps_frame = ttk.LabelFrame(self.master, text='Phase maps')
        self.import_map_button = ttk.Button(self.import_maps_frame,
                                            text='Import map', command=self.import_map_callback)
        self.clear_list_button = ttk.Button(self.import_maps_frame, text='Clear', command=self.clear_maps)

        self.maps_var = StringVar()
        self.maps_var.set('')
        if len(self.SLM.maps) > 0:
            self.maps = [m for m in self.SLM.maps]
        else:
            self.maps = ['Zeros']
        self.maps_options = OptionMenu(self.import_maps_frame, self.maps_var, *self.maps)

        self.maps_options.grid(column=0, row=0)
        self.import_map_button.grid(column=1, row=0)
        self.clear_list_button.grid(column=1, row=1)

        # ============================================================================================
        # Set up center(s) position
        # =============================================================================================
        # default mouse position for center is center of SLM
        self.mouse_coordinates = (int(self.SLM.width/2), int(self.SLM.height/2))
        self.center_position = [[int(self.SLM.width/2), int(self.SLM.height/2)]]
        self.plot_update()
        self.center_step = 1

        # =============================================================================================
        # Phase mask activation/de-activation
        # =============================================================================================
        self.active_frame = LabelFrame(self.master, text='Activate')
        self.active_var = StringVar()
        self.active_var.set('OFF')
        self.activation_button = Button(self.active_frame, textvariable=self.active_var,
                                        command=self.activate,  bg='firebrick2')
        self.activation_button.grid(column=0, row=0)
        self.active = False

        # ==========================================================================================
        # OPTIONS FRAME
        # ==========================================================================================

        self.notebook = ttk.Notebook(self.master)
        self.fqpm_frame = Frame(self.notebook)
        self.vortex_frame = Frame(self.notebook)
        self.multiple_frame = Frame(self.notebook)
        self.zernike_frame = Frame(self.notebook)
        self.notebook.add(self.fqpm_frame, text='FQ/EO')
        self.notebook.add(self.vortex_frame, text='Vortex')
        self.notebook.add(self.multiple_frame, text='Mutliple')
        self.notebook.add(self.zernike_frame, text='Zernike')
        self.notebook.grid()

        # ===========================================================================================
        # Star info in multiple star frame
        # ===========================================================================================

        self.stars_frame = ttk.LabelFrame(self.multiple_frame, text='Stars')
        self.star_1 = Label(self.stars_frame, text='Star 1')
        self.star_2 = Label(self.stars_frame, text='Star 2', state=DISABLED)
        self.star_3 = Label(self.stars_frame, text='Star 3', state=DISABLED)
        self.star_1.grid(column=0, row=1)
        self.star_2.grid(column=0, row=2)
        self.star_3.grid(column=0, row=3)

        I_lab = ttk.Label(self.stars_frame, text='Intensity', width=10)
        magn_lab = ttk.Label(self.stars_frame, text='Magnitude', width=10)
        l_lab = ttk.Label(self.stars_frame, text='Wavelength(nm)', width=10)
        F_lab = ttk.Label(self.stars_frame, text='F #', width=10)
        lD_lab = ttk.Label(self.stars_frame, text='l/D', width=10)
        phi_lab= ttk.Label(self.stars_frame, text='phi(pi)', width=10)
        C_lab = ttk.Label(self.stars_frame, text='Center(x,y)', width=10)
        magn_lab.grid(column=1, row=0)
        I_lab.grid(column=2, row=0)
        l_lab.grid(column=3, row=0)
        F_lab.grid(column=4, row=0)
        lD_lab.grid(column=5, row=0)
        phi_lab.grid(column=6, row=0)
        C_lab.grid(column=7, row=0)

        # 1st star -- always visible

        self.M1 = StringVar()
        self.M1.set('0')
        M1_entry = ttk.Entry(self.stars_frame, textvariable=self.M1, width=10)
        M1_entry.grid(column=1, row=1)
        self.I1_num = StringVar()
        self.I1_num.set('1')
        self.I1_entry = ttk.Entry(self.stars_frame, textvariable=self.I1_num, width=10)
        self.I1_entry.grid(column=2, row=1)
        self.l1_num = StringVar()
        self.l1_num.set('633')
        self.l1_entry = ttk.Entry(self.stars_frame, textvariable=self.l1_num, width=10)
        self.l1_entry.grid(column=3, row=1)
        self.F1_num = StringVar()
        self.F1_num.set('230')
        self.F1_entry = ttk.Entry(self.stars_frame, textvariable=self.F1_num, width=10)
        self.F1_entry.grid(column=4, row=1)
        self.starc1 = StringVar()
        self.starc1.set('%i,%i' % (int(self.SLM.width/2), int(self.SLM.height/2)))
        self.center1_lab = Entry(self.stars_frame, textvariable=self.starc1, width=10)
        self.center1_lab.grid(column=7, row=1)
        # star 2
        self.M2 = StringVar()
        self.M2.set('0')
        self.M2_entry = ttk.Entry(self.stars_frame, textvariable=self.M2,
                         width=10, state=DISABLED)
        self.M2_entry.grid(column=1, row=2)
        self.M2_entry.bind("<Return>", self.magnitude_to_intensity)
        self.I2_num = StringVar()
        self.I2_num.set('1')
        self.I2_entry = ttk.Entry(self.stars_frame, textvariable=self.I2_num,
                                  width=10, state=DISABLED)
        self.I2_entry.bind("<Return>", self.magnitude_to_intensity)
        self.I2_entry.grid(column=2, row=2)
        self.l2_num = StringVar()
        self.l2_num.set('633')
        self.l2_entry = ttk.Entry(self.stars_frame, textvariable=self.l2_num,
                                  width=10, state=DISABLED)
        self.l2_entry.grid(column=3, row=2)
        self.F2_num = StringVar()
        self.F2_num.set('230')
        self.F2_entry = ttk.Entry(self.stars_frame, textvariable=self.F2_num,
                                  width=10, state=DISABLED)
        self.F2_entry.grid(column=4, row=2)
        self.starc2 = StringVar()
        self.starc2.set('0,0')
        self.lD_star2 = StringVar()
        self.lD_star2.set('1')
        self.lD_star2_entry = Entry(self.stars_frame, textvariable=self.lD_star2,
                                    width=10, state=DISABLED)
        self.lD_star2_entry.grid(column=5, row=2)
        self.phi_star2 = StringVar()
        self.phi_star2.set('0')
        self.phi_star2_entry = Entry(self.stars_frame, textvariable=self.phi_star2,
                                     width=10, state=DISABLED)
        self.phi_star2_entry.grid(column=6, row=2)
        self.center2_lab = Entry(self.stars_frame, textvariable=self.starc2,
                                 width=10, state=DISABLED)
        self.center2_lab.grid(column=7, row=2)
        self.center2_lab.bind("<Return>", self.l_over_D_callback)
        # star 3
        self.M3 = StringVar()
        self.M3.set('0')
        self.M3_entry = ttk.Entry(self.stars_frame, textvariable=self.M3,
                         width=10, state=DISABLED)
        self.M3_entry.grid(column=1, row=3)
        self.M3_entry.bind("<Return>", self.magnitude_to_intensity)
        self.I3_num = StringVar()
        self.I3_num.set('1')
        self.I3_entry = ttk.Entry(self.stars_frame, textvariable=self.I3_num,
                                  width=10, state=DISABLED)
        self.I3_entry.grid(column=2, row=3)
        self.I3_entry.bind("<Return>", self.magnitude_to_intensity)
        self.l3_num = StringVar()
        self.l3_num.set('633')
        self.l3_entry = ttk.Entry(self.stars_frame, textvariable=self.l3_num,
                                  width=10, state=DISABLED)
        self.l3_entry.grid(column=3, row=3)
        self.F3_num = StringVar()
        self.F3_num.set('230')
        self.F3_entry = ttk.Entry(self.stars_frame, textvariable=self.F3_num,
                                  width=10, state=DISABLED)
        self.F3_entry.grid(column=4, row=3)
        self.starc3 = StringVar()
        self.starc3.set('0,0')
        self.lD_star3 = StringVar()
        self.lD_star3.set('1')
        self.lD_star3_entry = Entry(self.stars_frame, textvariable=self.lD_star3,
                                    width=10, state=DISABLED)
        self.lD_star3_entry.grid(column=5, row=3)
        self.phi_star3 = StringVar()
        self.phi_star3.set('0')
        self.phi_star3_entry = Entry(self.stars_frame, textvariable=self.phi_star3,
                                     width=10, state=DISABLED)
        self.phi_star3_entry.grid(column=6, row=3)
        self.center3_lab = Entry(self.stars_frame, textvariable=self.starc3,
                                 width=10, state=DISABLED)
        self.center3_lab.grid(column=7, row=3)
        self.center3_lab.bind("<Return>", self.l_over_D_callback)

        # ============================================================================================
        # FQPM and EOPM frame
        # ============================================================================================

        self.center1_lab_fqpm = Entry(self.fqpm_frame, textvariable=self.starc1)
        self.center1_lab_fqpm.grid(column=4, row=0)

        self.single_button = ttk.Button(self.fqpm_frame, text='Make map',
                                        command=lambda: self.make_map('single'))
        self.single_button.grid(column=0, row=0)
        map_types = ['FQPM', 'EOPM', 'FLAT']
        self.map_type_var = StringVar()
        self.map_type_var.set('FQPM')
        self.map_type_menu = OptionMenu(self.fqpm_frame, self.map_type_var, *map_types)
        self.map_type_menu.grid(row=0, column=2)

        # =========================================================================================================
        # CONTROL FRAME
        # =========================================================================================================
        self.control_frame = ttk.LabelFrame(self.master, text='Center Controls')
        self.cstep_var = StringVar()
        self.cstep_var.set('1')
        self.center_step_entry = Entry(self.control_frame, textvariable=self.cstep_var, justify='center')
        self.center_step_entry.bind("<Return>", self.set_center_step)
        self.center_control_up = ttk.Button(self.control_frame, text='^', command=lambda: self.center_move('up', 0))
        self.center_control_down = ttk.Button(self.control_frame, text='v', command=lambda: self.center_move('down',0))
        self.center_control_left = ttk.Button(self.control_frame, text='<', command=lambda: self.center_move('left',0))
        self.center_control_right = ttk.Button(self.control_frame, text='>', command=lambda: self.center_move('right',0))
        self.center_control_up.grid(column=1, row=0)
        self.center_control_down.grid(column=1, row=2)
        self.center_control_left.grid(column=0, row=1)
        self.center_control_right.grid(column=2, row=1)
        self.center_step_entry.grid(column=1, row=1)
        self.center_num = ['1']
        self.center_var = StringVar()
        self.center_var.set('1')

        # Set gray values
        self.val_1 = 0
        self.val_2 = 1
        self.grayval_frame = ttk.LabelFrame(self.fqpm_frame, text='Gray values')

        self.gray_1_val = StringVar()
        self.gray_1_val.set('0')
        self.gray_1_entry = Entry(self.grayval_frame, textvariable=self.gray_1_val, justify='center')
        self.gray_1_entry.bind("<Return>", self.arrow_return)
        self.gray_1_entry.bind("<Up>", self.arrow_return)
        self.gray_1_entry.bind("<Down>", self.arrow_return)
        self.gray_1_entry.bind("<Left>", self.arrow_return)
        self.gray_1_entry.bind("<Right>", self.arrow_return)
        self.gray_2_val = StringVar()
        self.gray_2_val.set('0')
        self.gray_2_entry = Entry(self.grayval_frame, textvariable=self.gray_2_val, justify='center')
        self.gray_2_entry.bind("<Return>", self.arrow_return)
        self.gray_2_entry.bind("<Up>", self.arrow_return)
        self.gray_2_entry.bind("<Down>", self.arrow_return)
        self.gray_2_entry.bind("<Left>", self.arrow_return)
        self.gray_2_entry.bind("<Right>", self.arrow_return)
        self.gray_1_lab = ttk.Label(self.grayval_frame, text='Gray-value 1')
        self.gray_2_lab = ttk.Label(self.grayval_frame, text='Gray-value 2')
        self.phase_1_val = StringVar()
        self.phase_1_val.set('Phase: %.3f rad'%self.menu.phase_curve[int(self.gray_1_val.get())])
        self.phase_2_val = StringVar()
        self.phase_2_val.set('Phase: %.3f rad'%self.menu.phase_curve[int(self.gray_2_val.get())])
        self.phase_1_lab = ttk.Label(self.grayval_frame, textvariable=self.phase_1_val)
        self.phase_2_lab = ttk.Label(self.grayval_frame, textvariable=self.phase_2_val)
        self.gray_1_lab.grid(column=0, row=0)
        self.gray_2_lab.grid(column=0, row=1)
        self.gray_1_entry.grid(column=1, row=0)
        self.gray_2_entry.grid(column=1, row=1)
        self.phase_1_lab.grid(column=2, row=0)
        self.phase_2_lab.grid(column=2, row=1)

        # ============================================================================================
        # ZERNIKE TAB
        # ============================================================================================
        # shift_per_pixel_lab = ttk.Label(self.fqpm_frame, text='Shift per px (horizontally)')
        # shift_per_pixel_lab.grid(column=0, row=3)
        # self.shift_per_pixel = DoubleVar()
        # self.shift_per_pixel.set(0.0)
        # shift_per_pixel_entry = Entry(self.fqpm_frame, textvariable=self.shift_per_pixel)
        # shift_per_pixel_entry.grid(column=1, row=3)

        defocus_coeff_lab = ttk.Label(self.zernike_frame, text='Defocus coeff:')
        defocus_coeff_lab.grid(column=0, row=0)
        self.defocus_coeff = DoubleVar()
        self.defocus_coeff.set(0)
        defocus_coeff_entry = Entry(self.zernike_frame, textvariable=self.defocus_coeff)
        defocus_coeff_entry.grid(column=1, row=0)

        astigm_coeff_lab = ttk.Label(self.zernike_frame, text='Astigmatism coeff:')
        astigm_coeff_lab.grid(column=2, row=0)
        self.astigm_coeff = DoubleVar()
        self.astigm_coeff.set(0)
        astigm_coeff_entry = Entry(self.zernike_frame, textvariable=self.astigm_coeff)
        astigm_coeff_entry.grid(column=3, row=0)

        zernike_range_lab = Label(self.zernike_frame, text='Phase shift of zernike')
        zernike_range_lab.grid(column=0, row=1, columnspan=2)
        self.zernike_min = DoubleVar()
        self.zernike_min.set(0)
        zernike_min_entry = Entry(self.zernike_frame, textvariable=self.zernike_min)
        zernike_min_entry.grid(column=2, row=1)
        self.zernike_max = DoubleVar()
        self.zernike_max.set(1)
        zernike_max_entry = Entry(self.zernike_frame, textvariable=self.zernike_max)
        zernike_max_entry.grid(column=3, row=1)

        apply_zernike = ttk.Button(self.zernike_frame, text='Apply', command=self.apply_zernike)
        apply_zernike.grid(column=0, row=2)
        self.Defocus = lambda r: np.sqrt(3)*(2*r**2)
        self.Astigm = lambda r, theta:np.sqrt(6)*r**2*np.sin(2*theta)

        xx, yy = np.meshgrid(np.arange(-self.SLM.width/2, self.SLM.width/2),
                             np.arange(-self.SLM.height/2, self.SLM.height/2))

        self.R, self.Theta = cart2pol(xx, yy)
        # ======================================================================================
        self.grayval_frame.grid(column=0, row=1, columnspan=5)
        self.control_frame.grid(column=0, row=2, columnspan=5)

        # ======================================================================================
        # Multiple sources
        # ======================================================================================
        # Pack star center vars together for easy access
        self.center_labels = [self.starc1, self.starc2, self.starc3]

        # make frame where a binary star map or triple star map can be created
        # binary phase map using airy pattern distribution for each star
        self.binary_frame = ttk.Frame(self.multiple_frame)
        self.binary_button = ttk.Button(self.binary_frame, text='Binary',
                                        command=lambda: self.make_map('binary'), state=DISABLED)
        self.binary_button.grid(column=1, row=1)
        self.checkbox_val = IntVar()
        binary_checkbox = Checkbutton(self.binary_frame, text='Save map', variable=self.checkbox_val)
        binary_checkbox.grid(column=3, row=1)

        self.tertiary_button = ttk.Button(self.binary_frame, text='Tertiary star',
                                        command=lambda: self.make_map('triple'), state=DISABLED)
        self.tertiary_button.grid(column=2, row=1)
        self.new_map_name = StringVar()
        self.new_map_name.set('Map name')
        self.new_map_name_entry = Entry(self.binary_frame, textvariable=self.new_map_name)
        self.new_map_name_entry_single = Entry(self.fqpm_frame, textvariable=self.new_map_name)
        self.new_map_name_entry_single.grid(column=3, row=0)
        self.new_map_name_entry.grid(column=0, row=1)
        self.save_filetypes = [('Windows Bitmap', '*.bmp'), ('Text File', '*.txt'), ('Fits File', '*.fits')]

        add_center_button = ttk.Button(self.binary_frame, text='Add', command=self.add_center)
        add_center_button.grid(column=0, row=0)
        self.centers_options = OptionMenu(self.binary_frame, self.center_var, *self.center_num)
        self.centers_options.grid(column=1, row=0)

        self.stars_frame.grid(column=0, row=0)
        self.binary_frame.grid(column=0, row=2)
        # =====================================================================================================
        # Vortex tab
        # =====================================================================================================
        self.make_vortex = ttk.Button(self.vortex_frame, text='Make vortex',
                                      command=lambda: self.make_map('vortex'))
        self.make_vortex.grid(column=0, row=0)
        # charge of the vortex
        charge_lab = ttk.Label(self.vortex_frame, text='charge')
        charge_lab.grid(column=2, row=1)
        self.charge = IntVar()
        self.charge.set(2)
        self.charge_entry = Entry(self.vortex_frame, textvariable=self.charge, width=10)
        self.charge_entry.bind("<Return>", self.charge_callback)
        self.charge_entry.grid(column=3, row=1)
        # coordinates entry
        coordinates_lab = ttk.Label(self.vortex_frame, text='Coordinates')
        coordinates_lab.grid(column=0, row=1)

        self.vortex_coordinates = StringVar()
        self.vortex_coordinates.set('%i, %i' % (int(self.SLM.width/2), int(self.SLM.height/2)))
        self.vortex_coordinates_entry = Entry(self.vortex_frame, textvariable=self.vortex_coordinates, width=10)
        self.vortex_coordinates_entry.grid(column=1, row=1)
        # label indicating gray values
        gray_lab = ttk.Label(self.vortex_frame, text='Gray values')
        gray_lab.grid(column=1, row=3, columnspan=2)
        # gray value for the 0 pi phase
        gray0_lab = ttk.Label(self.vortex_frame, text='0:', width=10)
        gray0_lab.grid(column=0, row=4)
        self.gray0 = IntVar()
        self.gray0.set(0)
        self.gray0_entry = Entry(self.vortex_frame, textvariable=self.gray0, width=10)
        self.gray0_entry.grid(column=1, row=4)
        # gray value for 2pi phase
        gray2pi_lab = ttk.Label(self.vortex_frame, text='2pi:', width=10)
        gray2pi_lab.grid(column=2, row=4)
        self.gray2pi = IntVar()
        self.gray2pi.set(0)
        self.gray2pi_entry = Entry(self.vortex_frame, textvariable=self.gray2pi, width=10)
        self.gray2pi_entry.grid(column=3, row=4)
        # button to change gray values of vortex on the fly
        self.gray_vortex_button = ttk.Button(self.vortex_frame, text='Change', command=self.vortex_change_grayvalues)
        self.gray_vortex_button.grid(column=4, row=4)
        # ======================================================================================================
        self.multiple_star_position = []
        # =============================================================================================================
        # Text frame
        # ========================================================================================================
        self.text_frame = ttk.Frame(self.master)
        scrollbar = Scrollbar(self.text_frame)
        scrollbar.grid(column=4, row=0)
        self.text = Text(self.text_frame, height=5, width=40, wrap='word', yscrollcommand=scrollbar.set)
        self.text.insert(INSERT, "Initializing SLM..\n")
        self.text.grid(column=0, row=0, columnspan=4)
        sys.stdout = StdoutRedirector(self.text)  # assign stdout to custom class

    def calculate_charge_gray(self, p):
        if not(0 <= self.gray0.get() <= 255 and 0 <= self.gray2pi.get() <= 255):
            print('invalid values')
            return
        if self.charge.get() % 2 != 0:
            print("Odd charge -> change to closest even number")
            self.charge.set(self.charge.get() + 1)
        #if 'vortex' not in self.SLM.maps.keys():
        #    return
        z = p**self.charge.get()
        z = (np.angle(z) + np.pi)/(2*np.pi)
        z = z*abs(self.gray2pi.get() - self.gray0.get()) + self.gray0.get()
        phase_map = np.zeros(self.SLM.dimensions, dtype=np.uint8)
        phase_map[:, :, 0] = z
        phase_map[:, :, 1] = z
        phase_map[:, :, 2] = z
        return phase_map

    def charge_callback(self, event):
        """
        Callback when charge of vortex phase mask is changed
        :return:
        """
        p = self.SLM.maps[self.maps_var.get()]['map']
        phase_map = self.calculate_charge_gray(p)
        self.image = phase_map
        print('Changed charge to %i' % self.charge.get())
        if self.active:
            self.SLM.draw(self.image)
        self.plot_update()
        return

    def make_vortex_callback(self):
        """
        Create single vortex mask at center denoted by star center
        :return:
        """
        try:
            c = self.vortex_coordinates.get().split(sep=',')
            xc = int(c[0])
            yc = int(c[1])
        except ValueError:
            print('Error with coordinates')
            return
        print('Calculating vortex with charge %i, gray %i-%i, coord %i,%i' %
              (self.charge.get(), self.gray0.get(), self.gray2pi.get(), xc, yc))

        p = self.SLM.Vortex_coronagraph(xc, yc)
        phase_map = self.calculate_charge_gray(p)

        name = "Vortex_coord:%i,%i" % (xc, yc)
        print('Finished, map-name %s' % name)
        self.SLM.maps[name] = {'data': phase_map}
        self.SLM.maps[name]['map'] = p
        self.SLM.maps[name]['center'] = [[xc, yc]]
        self.SLM.maps[name]['type'] = 'vortex'
        self.maps.append(name)
        self.refresh_optionmenu(self.maps_options, self.maps_var, self.maps)
        self.image = phase_map
        self.plot_update()
        # save map to bitmap if option is checked
        if self.checkbox_val.get():
            filename = filedialog.asksaveasfilename(parent=self.master, filetypes=self.save_filetypes,
                                                    title='Save map as..')
            filename += '.bmp'
            surf = pygame.surfarray.make_surface(phase_map)
            pygame.image.save(surf, filename)
        return

    def vortex_change_grayvalues(self):
        """
        Changes the values of the vortex by scaling them with new_range
        :return:
        """

        p = self.SLM.maps[self.maps_var.get()]['map']
        phase_map = self.calculate_charge_gray(p)
        self.image = phase_map
        print('Changed gray value range to %i-%i' % (self.gray0.get(), self.gray2pi.get()))
        if self.active:
            self.SLM.draw(self.image)
        self.plot_update()
        return

    def l_over_D_callback(self, event):
        """
        Transforms l/D and azimuthial information to pixel coordinates with respect to the first star
        x = n*l/D*cos(phi*pi)
        y = n*l/D*sin(phi*pi)
        :param which: which star
        :return:
        """
        if event.widget == self.center2_lab:
            x = int(float(self.lD_star2.get())*int(self.menu.lD.get())*np.cos(float(self.phi_star2.get())*np.pi-np.pi))
            y = int(float(self.lD_star2.get())*int(self.menu.lD.get())*np.sin(float(self.phi_star2.get())*np.pi-np.pi))
            self.multiple_star_position[0] = [float(self.lD_star2.get()), float(self.phi_star2.get())]
            x += self.center_position[0][0]
            y += self.center_position[0][1]
            self.center_position[1] = [x, y]
            self.multiple_star_position[0] = [x, y]
            self.starc2.set('%i,%i' % (x, y))
        elif event.widget == self.center3_lab:
            x = int(float(self.lD_star3.get())*int(self.menu.lD.get())*np.cos(float(self.phi_star3.get())*np.pi-np.pi))
            y = int(float(self.lD_star3.get())*int(self.menu.lD.get())*np.sin(float(self.phi_star3.get())*np.pi-np.pi))
            self.multiple_star_position[1] = [float(self.lD_star3.get()), float(self.phi_star3.get())]
            x += self.center_position[0][0]
            y += self.center_position[0][1]
            self.center_position[2] = [x, y]
            self.multiple_star_position[1] = [x, y]
            self.starc3.set('%i,%i' % (x, y))
        else:
            pass
        return

    def magnitude_to_intensity(self, event):
        """
        Transform magnitude difference between star and primary star into intensity difference
        :param which: which star
        :return:
        """
        if event.widget == self.M2_entry:
            m = float(self.M2.get())
            self.I2_num.set(str(10**(-m/2.5)))
        elif event.widget == self.M3_entry:
            m = float(self.M3.get())
            self.I2_num.set(str(10**(-m/2.5)))
        elif event.widget == self.I2_entry:
            I = float(self.I2_num.get())
            self.M2.set(-2.5*np.log10(I))
        elif event.widget == self.I3_entry:
            I = float(self.I3_num.get())
            self.M3.set(-2.5*np.log10(I))
        else:
            pass

        return

    def clear_maps(self):
        """
        Clears list of maps
        :return:
        """
        self.maps = []
        self.refresh_optionmenu(self.maps_options, self.maps_var, self.maps)
        return

    def make_map(self, which):
        """
        Make thread that starts to calculate phase map
        :param which:
        :return:
        """
        if which == 'single':
            self.map_thread = threading.Thread(target=self.single_mask, daemon=True)
            self.map_thread.start()
        elif which == 'binary':
            self.map_thread = threading.Thread(target=self.binary_mask, daemon=True)
            self.map_thread.start()
        elif which == 'triple':
            self.map_thread = threading.Thread(target=self.triple_mask, daemon=True)
            self.map_thread.start()
        elif which == 'vortex':
            self.map_thread = threading.Thread(target=self.make_vortex_callback, daemon=True)
            self.map_thread.start()
        else:
            pass

        print('Thread started')
        return

    def refresh_optionmenu(self, menu, var, options):
        """
        Refreshes option menu
        :param menu: handle to optionmenu widget
        :param var: handle to variable of menu
        :param options: options to insert
        :return:
        """
        var.set('')
        menu['menu'].delete(0, 'end')
        if len(options) == 0:
            menu['menu'].add_command(label='', command=_setit(var, ''))
            return
        for option in options:
            menu['menu'].add_command(label=option, command=_setit(var, option))
        return

    def add_center(self):
        """
        Add new center which represents a new star or other object on the phase map
        Gets center coordinates from right click mouse position on figure.
        :return:
        """
        if len(self.center_num) > 2:  # up to a total of 3 objects can be defined
            return
        num = int(self.center_num[-1]) + 1
        self.center_num.append(str(num))
        self.center_position.append([0, 0])
        self.multiple_star_position.append([1, 0])  # [l/D, phi]
        self.refresh_optionmenu(self.centers_options, self.center_var, self.center_num)
        if num == 2:
            self.M2_entry.config(state=NORMAL)
            self.I2_entry.config(state=NORMAL)
            self.l2_entry.config(state=NORMAL)
            self.F2_entry.config(state=NORMAL)
            self.lD_star2_entry.config(state=NORMAL)
            self.phi_star2_entry.config(state=NORMAL)
            self.star_2.config(state=NORMAL)
            self.binary_button.config(state=NORMAL)
            self.center2_lab.config(state=NORMAL)
            self.multiple_star_position.append([0, 0])
        else:
            self.M3_entry.config(state=NORMAL)
            self.I3_entry.config(state=NORMAL)
            self.l3_entry.config(state=NORMAL)
            self.F3_entry.config(state=NORMAL)
            self.lD_star3_entry.config(state=NORMAL)
            self.phi_star3_entry.config(state=NORMAL)
            self.tertiary_button.config(state=NORMAL)
            self.star_3.config(state=NORMAL)
            self.center3_lab.config(state=NORMAL)
            self.multiple_star_position.append([0, 0])
        return

    def arrow_return(self, event):
        """
        Changes grayvalue with arrow key hits
        :param event:
        :return:
        """
        which = event.widget
        what = event.keycode
        if what == 37 or what == 40:
            what = -1
        elif what == 38 or what == 39:
            what = 1
        elif what == 13:
            what = 0
        else:
            return
        if which == self.gray_1_entry:
            try:
                val_old = self.val_1
                val = int(self.gray_1_val.get())
                val += what
                if val > 255:
                    val = 255
                if val < 0:
                    val = 0
                self.gray_1_val.set(str(val))
                self.phase_1_val.set('Phase: %.3f rad'%self.menu.phase_curve[val])
                self.val_1 = val
                self.set_value(val_old, val)
            except ValueError:
                return
        elif which == self.gray_2_entry:
            try:
                val_old = self.val_2
                val = int(self.gray_2_val.get())
                val += what
                if val > 255:
                    val = 255
                if val < 0:
                    val = 0
                self.gray_2_val.set(str(val))
                self.phase_2_val.set('Phase: %.3f rad'%self.menu.phase_curve[val])
                self.val_2 = val
                self.set_value(val_old, val)
            except ValueError:
                return
        else:
            return

    def set_value(self, val_old, val):
        """
        Find all pixels with value val and replace them with the new one
        :param val_old: old value to replace
        :param val: value to replace with
        :return:
        """
        self.image[self.image == val_old] = val
        if self.active:
            self.SLM.draw(self.image)
        self.plot_update()
        return

    def activate(self):
        """
        Activate and deactivate SLM
        :return:
        """
        if self.active:
            self.active = False
            self.activation_button.config(bg='firebrick2')
            self.active_var.set('OFF')
            self.send_map('OFF')
        else:
            self.active = True
            self.activation_button.config(bg='PaleGreen2')
            self.active_var.set('ON')
            self.send_map('ON')
        return

    def get_grayvals(self, image):
        """
        Get the values from the phase mask
        :param image: applied mask
        :return:
        """
        vals = np.unique(image)
        self.val_1 = vals.min()
        self.val_2 = vals.max()
        self.gray_1_val.set(str(self.val_1))
        self.gray_2_val.set(str(self.val_2))
        return

    def send_map(self, status):
        """
        Send map to SLM
        :param status: Phase map for ON and zero map for OFF
        :return:
        """
        if status == 'ON':
            map_array = self.SLM.maps[self.maps_var.get()]['data']  # should get the matrix of the chosen map
            self.get_grayvals(map_array)
            self.image = map_array
            self.SLM.draw(map_array)
            self.plot_update()
        elif status == 'OFF':
            self.SLM.maps[self.maps_var.get()]['data'] = self.image  # save current state of map to dictionary
            self.image = self.off_image
            self.SLM.draw(self.off_image)
            self.plot_update()

    def plot_update(self):
        self.im.set_data(self.image[:, :, 0].T)
        self.fig.canvas.draw()
        return

    def import_map_callback(self):
        """
        Import new map from file. Accepted extensions are bmp, txt, fits
        :return:
        """
        mfile = filedialog.askopenfilename()
        try:
            mname = self.SLM.import_phase_map(mfile)
            mname = os.path.basename(mname)
            self.maps.append(mname)
            self.refresh_optionmenu(self.maps_options, self.maps_var, self.maps)
        except Exception as e:
            print(e)
            return

        return

    def click_callback(self, event):
        _x = event.x
        _y = event.y
        inv = self.ax.transData.inverted()
        data_pos = inv.transform((_x, _y))
        data_pos = tuple([int(e) for e in data_pos])
        if event.button == 1:
            self.center_move('mouse', data_pos)
        elif event.button == 3:
            self.mouse_coordinates = data_pos
        else:
            pass
        return

    def center_move(self, d, pos):
        """
        Move center of phase map
        :param d: direction to move
        :return:
        """
        which = int(self.center_var.get())-1  # which center is currently active
        if d == 'up':
            self.center_position[which][1] -= self.center_step
            self.image = np.roll(self.image, shift=-self.center_step, axis=1)
        elif d == 'down':
            self.center_position[which][1] += self.center_step
            self.image = np.roll(self.image, shift=self.center_step, axis=1)
        elif d == 'left':
            self.center_position[which][0] -= self.center_step
            self.image = np.roll(self.image, shift=-self.center_step, axis=0)
        elif d == 'right':
            self.center_position[which][0] += self.center_step
            self.image = np.roll(self.image, shift=self.center_step, axis=0)
        else:
            self.center_position[which] = list(pos)
        # update plot and SLM
        if self.active:
            self.SLM.draw(self.image)
        self.plot_update()
        # update label showing center
        s = '%i,%i'%(self.center_position[which][0], self.center_position[which][1])
        self.center_labels[which].set(s)
        return

    def set_center_step(self, event):
        """
        Callback for setting the center move step size
        """
        try:
            val = int(self.cstep_var.get())
            if val > 384:
                raise ValueError
            self.center_step = val
        except ValueError:
            self.cstep_var.set(str(self.center_step))
        return

    def binary_mask(self):
        """
        Create binary mask
        :return:
        """

        c1 = self.starc1.get().split(sep=',')
        c1 = (int(c1[0]), int(c1[1]))
        c2 = self.starc2.get().split(sep=',')
        c2 = (int(c2[0]), int(c2[1]))
        try:
            I1, l1, F1 = float(self.I1_num.get()), float(self.l1_num.get()), float(self.F1_num.get())
            I2, l2, F2 = float(self.I2_num.get()), float(self.l2_num.get()), float(self.F2_num.get())
            val1 = self.menu.phase_curve[int(self.gray_1_val.get())]
            val2 = self.menu.phase_curve[int(self.gray_2_val.get())]
            print('Binary map with values :%f, %f'%(val1, val2))
        except ValueError:
            print('ValueError')
            return
        self.f = lambda x, y: self.SLM.pixel_value(x, y, c1, c2, I1, I2, val1, val2, F1, F2, l1, l2,
                                                   mask=self.map_type_var.get())
        print('Calculating binary %s' % self.map_type_var.get())
        p = np.zeros(np.SLM.size)
        print('Running binary weight-values calculation..')
        for (x, y), val in np.ndenumerate(p):
            p[x, y] = self.f(x, y)

        try:
            print("Running rad to gray conversion..")
            p = self.rad_to_gray(p)  # convert rad values to the corresponding gray values
        except Exception as e:
            print(e)
            return

        phase_map = np.zeros(self.SLM.dimensions, dtype=np.uint8)  # dimensions are (width,height, 3)
        phase_map[:, :, 0] = p
        phase_map[:, :, 1] = p
        phase_map[:, :, 2] = p
        name = self.new_map_name.get()
        self.SLM.maps[name] = {'data': phase_map}
        self.SLM.maps[name]['center'] = [[c1[0], c1[1]], [c2[0], c2[1]]]
        self.SLM.maps[name]['star info'] = [[I1, l1, F1], [I2, l2, F2]]
        self.maps.append(name)
        self.refresh_optionmenu(self.maps_options, self.maps_var, self.maps)
        self.image = phase_map
        self.plot_update()
        # save map to bitmap if option is checked
        if self.checkbox_val.get():
            self.menu.save_fits(name=name)
            print('File saved')
        return

    def triple_mask(self):
        """
        Create binary mask
        :return:
        """

        c1 = self.starc1.get().split(sep=',')
        c1 = (int(c1[0]), int(c1[1]))
        c2 = self.starc2.get().split(sep=',')
        c2 = (int(c2[0]), int(c2[1]))
        c3 = self.starc3.get().split(sep=',')
        c3 = (int(c3[0]), int(c3[1]))

        try:
            I1, l1, F1 = int(self.I1_num.get()), int(self.l1_num.get()), int(self.F1_num.get())
            I2, l2, F2 = int(self.I2_num.get()), int(self.l2_num.get()), int(self.F2_num.get())
            I3, l3, F3 = int(self.I3_num.get()), int(self.l3_num.get()), int(self.F3_num.get())
            val1 = int(self.gray_1_val.get())
            val2 = int(self.gray_2_val.get())
        except ValueError:
            print('Error')
            return
        #self.f = lambda x, y: self.SLM.pixel_value_triple(x, y, c1, c2, I1, I2, val1, val2, F1, F2, l1, l2)
        p = np.zeros(self.SLM.size, dtype=np.uint8)

        for (x, y), val in np.ndenumerate(p):
            p[x, y] = self.f(x, y)

        phase_map = np.zeros(self.SLM.dimensions, dtype=np.uint8)
        phase_map[:, :, 0] = p
        phase_map[:, :, 1] = p
        phase_map[:, :, 2] = p
        name = self.new_map_name.get()
        self.SLM.maps[name] = {'data': phase_map}
        self.SLM.maps[name]['center'] = [[c1[0], c1[1]], [c2[0], c2[1]], [c3[0], c3[1]]]
        self.SLM.maps[name]['star info'] = [[I1, l1, F1], [I2, l2, F2], [I3, l3, F3]]
        self.maps.append(name)
        self.refresh_optionmenu(self.maps_options, self.maps_var, self.maps)
        self.image = phase_map
        self.plot_update()
        # save map to bitmap if option is checked
        if self.checkbox_val.get():
            self.menu.save_fits(name=name)
            print('File saved')
        return

    def single_mask(self):
        """
        Create single star mask at center denoted by star center
        :return:
        """
        try:
            which = int(self.center_var.get())-1
            c = self.center_labels[which].get().split(sep=',')
            I1, l1, F1 = int(self.I1_num.get()), int(self.l1_num.get()), int(self.F1_num.get())
            xc = int(c[0])
            yc = int(c[1])
            val1 = int(self.gray_1_val.get())
            val2 = int(self.gray_2_val.get())
        except ValueError:
            print('Error')
            return
        p = np.zeros(self.SLM.size)
        p_raw = np.zeros(self.SLM.size)
        print('Calculating %s with gray values %i, %i at coord %i,%i' %
              (self.map_type_var.get(), val1, val2, xc, yc))
        if self.map_type_var.get() == 'FQPM':
            p[xc:, yc:] = val2
            p[:xc, :yc] = val2
            p[:xc, yc:] = val1
            p[xc:, :yc] = val1
            p_raw[xc:, yc:] = float(val2/255)*2*np.pi
            p_raw[:xc, :yc] = float(val2/255)*2*np.pi
            p_raw[:xc, yc:] = float(val1/255)*2*np.pi
            p_raw[xc:, :yc] = float(val1/255)*2*np.pi
        elif self.map_type_var.get() == 'FLAT':
            pass
        elif self.map_type_var.get() == 'EOPM':
            for (x, y), v in np.ndenumerate(p):
                p[x, y] = self.SLM.eight_octants(x, y, (xc, yc), val1, val2)
                p_raw = p
        else:
            # would be nice to have some feedback in the GUI at some point
            return
        phase_map = np.zeros(self.SLM.dimensions, dtype=np.uint8)
        phase_map[:, :, 0] = p
        phase_map[:, :, 1] = p
        phase_map[:, :, 2] = p
        name = self.new_map_name.get()
        self.SLM.maps[name] = {'data': phase_map}
        self.SLM.maps[name]['center'] = [[xc, yc]]
        self.SLM.maps[name]['star info'] = [I1, l1, F1]
        self.SLM.maps[name]['map'] = p_raw
        self.maps.append(name)
        self.refresh_optionmenu(self.maps_options, self.maps_var, self.maps)
        print('Finished, mask-name: %s' % name)
        self.image = phase_map
        self.plot_update()
        # save map to bitmap if option is checked
        if self.checkbox_val.get():
            filename = filedialog.asksaveasfilename(parent=self.master, filetypes=self.save_filetypes,
                                                    title='Save map as..')
            filename += '.bmp'
            surf = pygame.surfarray.make_surface(phase_map)
            pygame.image.save(surf, filename)
            print('File saved')
        return

    def zernike_send(self):
        """
        Changes the values of the vortex by scaling them with new_range
        :return:
        """

        p = self.SLM.maps[self.maps_var.get()]['data']

        self.image = p
        if self.active:
            self.SLM.draw(self.image)
        self.plot_update()
        return

    def apply_zernike(self):
        zernike = self.defocus_coeff.get()*self.Defocus(self.R/1080) + \
                  self.astigm_coeff.get()*self.Astigm(self.R/1080, self.Theta)
        magnitude = (self.zernike_max.get()-self.zernike_min.get())
        p = self.SLM.maps[self.maps_var.get()]['map']
        p = np.exp(1j*p)
        calib = np.angle(np.exp(1j*zernike.T*magnitude))
        m = (np.angle(p/np.exp(1j*calib)) + np.pi)/(2*np.pi)
        m *= 255
        phase_map = np.zeros(self.SLM.dimensions, dtype=np.uint8)
        phase_map[:, :, 0] = m
        phase_map[:, :, 1] = m
        phase_map[:, :, 2] = m
        self.SLM.maps[self.maps_var.get()]['data'] = phase_map
        self.zernike_send()
        return

    def rad_to_gray(self, p):
        """
        Transforms the map p which contains values in radians to the nearest fitting gray values
        :param p: map in radians
        :return: transformed map
        """
        # first find the unique values of p
        p = np.round(self.menu.rad_2_gray(p))
        p[p < 0] = 0  # correct possible negative values of fit
        p[p > 255] = 255  # if value gets over 255
        return p.astype(np.uint8)


if __name__ == '__main__':
    def on_closing(master):
        master.quit()
        master.destroy()
        return

    root = Tk()
    window = SLMViewer(root)
    window.data_plot.get_tk_widget().grid(column=0, row=0, columnspan=4, rowspan=5)
    window.import_maps_frame.grid(column=4, row=5)
    window.text_frame.grid(column=4, row=3)
    window.control_frame.grid(column=4, row=2)
    #window.grayval_frame.grid(column=4, row=2)
    window.active_frame.grid(column=4, row=0)
    #window.stars_frame.grid(column=4, row=3)
    #window.binary_frame.grid(column=4, row=4)
    #window.rotate_frame.grid(column=0, row=5)
    window.notebook.grid(column=4, row=1)
    root.protocol("WM_DELETE_WINDOW", lambda: on_closing(root))  # make sure window close properly
    root.mainloop()
