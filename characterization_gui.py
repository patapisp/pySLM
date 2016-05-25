__author__ = 'Chronis'
from definitions import SLM
import numpy as np
import pygame
from QtGUI.core.pco_definitions import PixelFly
from astropy.io import fits
import os, sys, time
from PyQt4 import QtCore, QtGui
import pyqtgraph as pg
from threading import Thread


class SLMInspector(QtGui.QWidget):
    """
    GUI that enables characterization of SLM. It provides a mutation of SLM and camera attributes that will facilitate
    obtaining characterization data/curves fot the SLM using the PixelFly camera.
    """
    def __init__(self, main_window, parent=None):
        #QtGui.QWidget.__init__(self, parent)
        super(SLMInspector, self).__init__()
        self.path = os.path.dirname(os.path.realpath("__file__"))
        try:
            self.camera = PixelFly(dllpath=self.path)
        except:
            print('dll SC2_Cam.dll not found')

        try:
            self.SLM = SLM()
        except:
            print('No second screen connected')

        self.camera_connected = False
        self.camera_alive = False

        self.slm_connected = False
        self.stack = []
        self.stack_names = []
        self.browse_index = 0
        self.central_widget = QtGui.QWidget(main_window)

        # Grid layout to place all widgets
        self.widget_layout = QtGui.QGridLayout()
        # layout to put images on
        self.gw_slm = pg.GraphicsLayoutWidget()
        # make margins around image items zero
        #self.gw_slm.ci.layout.setContentsMargins(0,0,0,0)

        main_window.setCentralWidget(self.central_widget)
        # connect camera button
        self.connect_camera_btn = QtGui.QPushButton('Connect Camera')
        self.connect_camera_btn.setStyleSheet("background-color: darkCyan")
        QtCore.QObject.connect(self.connect_camera_btn, QtCore.SIGNAL('clicked()'), self.connect_camera)
        self.widget_layout.addWidget(self.connect_camera_btn, 0, 0)
        self.vb = pg.ViewBox()
        self.gw_slm.addItem(self.vb)
        self.vb.setAspectLocked(lock=True, ratio=1)
        # invert Y axis -> PyQt <-> Numpy arrays convention
        self.vb.invertY()
        # Image Item is the image displaying item. Has a lot of options and the user can zoom in/out by pressing the
        # right mouse button and moving the mouse up/down. Furthermore by going over the image with the mouse will
        # indicate the coordinates and value.
        self.image = pg.ImageItem()
        self.image_text = pg.TextItem(text="Stack frame :", color=(200, 0, 0))
        """
        file = "binary_map.bmp"
        p = pygame.image.load(file)
        p = pygame.surfarray.array3d(p)
        """
        p = np.zeros((1024, 768, 3), dtype=np.uint8)
        self.image.setImage(p)

        self.vb.addItem(self.image)
        self.vb.addItem(self.image_text)
        self.widget_layout.addWidget(self.gw_slm, 1, 5, 5, 5)

        self.make_gray_stack_btn = QtGui.QPushButton('Make gray stack')
        gray_value_steps_lab = QtGui.QLabel('Step size:')
        self.gray_value_steps = QtGui.QLineEdit('10')

        self.import_stack_btn = QtGui.QPushButton('Import stack..')
        self.runs_lab = QtGui.QLineEdit('1')
        runs_label = QtGui.QLabel("Number of runs")
        self.sweep_stack_btn = QtGui.QPushButton('Sweep stack')
        number_of_frames_lab = QtGui.QLabel('Number of frames:')
        self.number_of_frames = QtGui.QLineEdit('10')
        self.measure_btn = QtGui.QPushButton('Measure')
        self.left_btn = QtGui.QPushButton("<-")
        self.right_btn = QtGui.QPushButton("->")
        QtCore.QObject.connect(self.measure_btn, QtCore.SIGNAL('clicked()'), self.measure_callback)
        QtCore.QObject.connect(self.right_btn, QtCore.SIGNAL('clicked()'), self.browse_stack)
        QtCore.QObject.connect(self.left_btn, QtCore.SIGNAL('clicked()'), self.browse_stack)
        QtCore.QObject.connect(self.make_gray_stack_btn, QtCore.SIGNAL('clicked()'), self.make_gray_stack)
        QtCore.QObject.connect(self.import_stack_btn, QtCore.SIGNAL('clicked()'), self.import_stack)
        QtCore.QObject.connect(self.sweep_stack_btn, QtCore.SIGNAL('clicked()'), self.sweep_return)
        self.widget_layout.addWidget(self.make_gray_stack_btn, 2, 0)
        self.widget_layout.addWidget(gray_value_steps_lab, 3, 0)
        self.widget_layout.addWidget(self.gray_value_steps, 3, 1)
        self.widget_layout.addWidget(self.import_stack_btn, 2, 1)
        self.widget_layout.addWidget(self.sweep_stack_btn, 4, 0)
        self.widget_layout.addWidget(self.measure_btn, 4, 1)
        self.widget_layout.addWidget(number_of_frames_lab, 5, 0)
        self.widget_layout.addWidget(self.runs_lab, 6, 1)
        self.widget_layout.addWidget(runs_label, 6, 0)
        self.widget_layout.addWidget(self.number_of_frames, 5, 1)
        self.widget_layout.addWidget(self.left_btn, 0, 6)
        self.widget_layout.addWidget(self.right_btn, 0, 7)
        self.central_widget.setLayout(self.widget_layout)
        self.central_widget.setFocusPolicy(QtCore.Qt.StrongFocus)
        self.central_widget.setFocus()

        self.save_dir = self.path
        self.weight = np.ones((1024, 768))
        self.weight[:, 384:] = 0

    def browse_stack(self):
        """
        By pressing arrow keys one can go through the stack image and view them in the plot
        :return:
        """

        if len(self.stack) == 0:
            print('Zero length')
            return
        if self.sender() == self.left_btn:
            if self.browse_index > 0:
                self.browse_index -= 1
            else:
                self.browse_index = 0
        elif self.sender() == self.right_btn:
            if self.browse_index < len(self.stack)-1:
                self.browse_index += 1
            else:
                self.browse_index = len(self.stack)-1
        mask = np.zeros((1024, 768, 3), dtype=np.uint8)
        mask[:,:,0]=self.stack[self.browse_index]*self.weight
        mask[:,:,1]=self.stack[self.browse_index]*self.weight
        mask[:,:,2]=self.stack[self.browse_index]*self.weight
        self.image.setImage(mask)
        self.image_text.setText('Stack frame : %s' % self.stack_names[self.browse_index], color=(200, 0, 0))
        self.SLM.draw(mask)
        return


    def connect_camera(self):
        """
        Connect to camera. If camera connection returns error report it and
        set connected status to False.
        :return:
        """
        if self.camera_connected:
            err = self.camera.close_camera()
            self.camera_connected = False
            self.connect_camera_btn.setText('CONNECT')
            self.connect_camera_btn.setStyleSheet("background-color: darkCyan")
            #self.connection_status.setText('Disconnected')
        else:
            err = self.camera.open_camera()
            if not err:
                #self.connection_status.setText('Error with connection')
                print('Error with connection')
                return
            self.camera_connected = True
            self.connect_camera_btn.setText('DISCONNECT')
            self.connect_camera_btn.setStyleSheet("background-color: green")
            #self.connection_status.setText('Connected')
        return

    def make_gray_stack(self):
        """
        Makes a stack of phase maps to go through with the whole map at one value for n
        values between 0-255
        n: number of gray values to use (max=255, min=2)
        :return: returns numpy array of n 1024x768x3 matrices
        """

        step = int(self.gray_value_steps.text())
        if step == 0 or step >= 255:
            return
        values = np.arange(0, 255, step)
        values = np.append(values, 255)
        self.stack = []
        self.stack_names = []
        for val in values:
            #self.stack.append(np.ones((1024, 768), dtype=np.uint8)*weight*val)
            self.stack_names.append('GrayValue_%03i' % (val,))
        self.stack = values
        self.browse_index = 0
        try:
            mask = np.ones((1024, 768, 3), dtype=np.uint8)
            mask[:,:,0]=self.stack[0]*self.weight
            mask[:,:,1]=self.stack[0]*self.weight
            mask[:,:,2]=self.stack[0]*self.weight
            self.image.setImage(mask)
        except IndexError:
            pass
        self.image_text.setText("Stack frame: %s" % self.stack_names[0], color=(200, 0, 0))
        return

    def import_stack(self):
        """
        Imports stack of frames from chosen folder. The functions expects the folder to contain only the phase maps in
        .fits format with a 1024x768 resolution. Furthermore, the files should be numbered with the _XXXX.fits
        convention.
        In the end each file to be added to the stack should look like /phase_map_name_0001.fits
        :return:
        """
        filename = QtGui.QFileDialog.getOpenFileName(self, 'Import stack..')

        name, extension = os.path.splitext(filename)
        if extension != '.fits':
            print('Error: Unrecognized extension %s' % extension)
            return
        path = os.path.dirname(filename)+'\\'
        self.stack = []
        self.stack_names = []
        files = [f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))]
        for f in files:
            save_name = os.path.splitext(f)[0]
            hdu = fits.open(os.path.join(path, f))
            p = hdu[0].data
            m = np.zeros((1024, 768, 3), dtype=np.uint8)
            m[:, :, 0] = p
            m[:, :, 1] = p
            m[:, :, 2] = p
            self.stack.append(m)
            self.stack_names.append(save_name)
            del hdu
        self.browse_index = 0
        self.image.setImage(self.stack[0])
        self.image_text.setText("Stack frame: %s" % self.stack_names[0], color=(200, 0, 0))
        return

    def sweep_return(self):
        if not self.camera_connected:
            return

        try:
            self.camera_alive = True
            self.camera.arm_camera()
            self.camera.start_recording()
            self.camera.allocate_buffer()
            self.camera._prepare_to_record_to_memory()
            # make thread to run image loop
            self.image_thread = Thread(target=self.camera.record_to_memory_2)
            self.image_thread.setDaemon(True)
            self.image_thread.start()
            self.filename = QtGui.QFileDialog.getSaveFileName(self, 'Save as..')
            self.run_index = 0
            self.runs = int(self.runs_lab.text())
            self.dirname = os.path.dirname(self.filename)
            d =self.dirname+ "\\run%i\\"%self.run_index
            os.makedirs(d)
            self.filename = os.path.join(d, os.path.basename(self.filename))

            self.meas_index = 0
            QtCore.QTimer.singleShot(500, self.timer_return)
        except:
            pass
        return

    def timer_return(self):
        if self.run_index < self.runs:
            num_of_frames = 10
            try:
                num_of_frames = int(self.number_of_frames.text())
            except ValueError:
                return
            if self.meas_index < len(self.stack):
                m  = self.stack[self.meas_index]
                name = self.stack_names[self.meas_index]
                mask = np.zeros((1024, 768, 3), dtype=np.uint8)
                mask[:,:,0]=m*self.weight
                mask[:,:,1]=m*self.weight
                mask[:,:,2]=m*self.weight
                self.SLM.draw(mask)
                time.sleep(1)
                #self.image.setImage(mask)
                name = self.filename + '_' + name
                print('Measuring value:%s'%name)
                self.record_and_save(name, num_of_frames)
                self.meas_index += 1
                QtCore.QTimer.singleShot(500, self.timer_return)
            else:
                self.run_index += 1
                self.meas_index = 0
                d =self.dirname+ "\\run%i\\"%self.run_index
                os.makedirs(d)
                self.filename = os.path.join(d, os.path.basename(self.filename))
                QtCore.QTimer.singleShot(500, self.timer_return)
        else:
            self.camera.live = False
            time.sleep(1)
            self.image_thread.join()
            self.camera.disarm_camera()
            self.camera_alive = False
            #self.sweep_thread.join()
            print('Finished sweep')





    def sweep_stack(self):
        """
        Sweep through currently active stack and take frames for each phase map on the stack. The number of frames if
        defined below the button.
        :return:
        """
        if not self.camera_connected:
            return
        weight = np.ones((1024, 768))
        weight[:, 384:] = 0
        try:
            self.camera_alive = True
            self.camera.arm_camera()
            self.camera.start_recording()
            self.camera.allocate_buffer()
            self.camera._prepare_to_record_to_memory()
            # make thread to run image loop
            self.image_thread = Thread(target=self.camera.record_to_memory_2)
            self.image_thread.setDaemon(True)
            self.image_thread.start()
            # ask for file name
            filename = QtGui.QFileDialog.getSaveFileName(self, 'Save as..')
            num_of_frames = 10
            try:
                num_of_frames = int(self.number_of_frames.text())
            except ValueError:
                return
            for m, name in zip(self.stack, self.stack_names):
                mask = np.zeros((1024, 768, 3), dtype=np.uint8)
                mask[:,:,0]=m*self.weight
                mask[:,:,1]=m*self.weight
                mask[:,:,2]=m*self.weight
                self.SLM.draw(mask)
                time.sleep(1)
                #self.image.setImage(mask)
                name = filename + '_' + name
                print('Measuring value:%s'%name)
                self.record_and_save(name, num_of_frames)
        # except:
        #     print('Error')
        except MemoryError:
            pass
        finally:
            self.camera.live = False
            time.sleep(1)
            self.image_thread.join()
            self.camera.disarm_camera()
            self.camera_alive = False
            #self.sweep_thread.join()
            print('Finished sweep')
        return

    def record_and_save(self, filename, frames):
        """
        Takes frames number of frames and stores them in .fits format under the name filename with a postfix
        filename_XXXX.fits
        :param filename: name of the file and path to save data
        :param frames: number of frames to save
        :return:
        """
        # check if camera is alive
        if not self.camera_alive:
            print('Camera not armed')
            return
        if not(isinstance(frames, int)):
            print('Not int number of frames')
            return
        d=0
        for j in range(2):
            d = self.get_queue_save(5)
            # flush memory so that correct frames are saved
        del d

        hdu = fits.PrimaryHDU()
        for i in range(frames):
            file = filename + "_%04d"%(i,)+'.fits'
            hdu.data = self.camera.q.get()
            # other header details will come in here
            hdu.writeto(file)
        return

    def measure_callback(self):
        if not self.camera_connected:
            return
        weight = np.ones((1024, 768))
        weight[:, 384:] = 0
        try:
            self.camera_alive = True
            self.camera.arm_camera()
            self.camera.start_recording()
            self.camera.allocate_buffer()
            self.camera._prepare_to_record_to_memory()
            # make thread to run image loop
            self.image_thread = Thread(target=self.camera.record_to_memory_2)
            self.image_thread.setDaemon(True)
            self.image_thread.start()
            # ask for file name
            filename = QtGui.QFileDialog.getSaveFileName(self, 'Save as..', self.save_dir)
            self.save_dir = os.path.dirname(filename)
            file = os.path.basename(filename)
            path = self.save_dir + "\\GV%s"%(self.stack_names[self.browse_index].split(sep='_')[-1])
            os.makedirs(path)
            num_of_frames = 10
            try:
                num_of_frames = int(self.number_of_frames.text())
            except ValueError:
                return
            for m, name in zip([self.stack[0], self.stack[self.browse_index]], [self.stack_names[0], self.stack_names[self.browse_index]]):
                mask = np.zeros((1024, 768, 3), dtype=np.uint8)
                mask[:,:,0]=m*self.weight
                mask[:,:,1]=m*self.weight
                mask[:,:,2]=m*self.weight
                self.SLM.draw(mask)
                time.sleep(1)
                #self.image.setImage(mask)

                name = path + "\\" + file + '_' + name
                print('Measuring value:%s'%name)
                self.record_and_save(name, num_of_frames)
        # except:
        #     print('Error')
        except MemoryError:
            pass
        finally:
            self.camera.live = False
            time.sleep(1)
            self.image_thread.join()
            self.camera.disarm_camera()
            self.camera_alive = False
            print('Finished measuring value')
        return


    def get_queue_save(self, frames):
        """
        Takes frames number of frames from the image queue
        :param frames: number of frames to get
        :return:
        """
        data = []
        i = 0
        while i < frames:
            try:
                data.append(self.camera.q.get())
                i += 1
            except:
                pass

        return data


if __name__ == '__main__':

    app = QtGui.QApplication(sys.argv)
    window = QtGui.QMainWindow()
    window.setWindowTitle('SLM Inspector')

    gui = SLMInspector(window)

    window.show()
    sys.exit(app.exec_())




