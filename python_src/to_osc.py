import tensorflow as tf

from pythonosc import udp_client
from pythonosc import dispatcher
from pythonosc import osc_server

import struct
import win32clipboard
import win32con
import win32gui
import win32process
import win32api
import os
import ctypes

IP = "127.0.0.1"
PORT_SEND = 5005
PORT_RECEIVE = 5006
MY_WINDOW_NAME = "Command Prompt"

def to_cf_hdrop(file_names):
    """
    Convert a list of files into CF_HDROP format,
    which is used to copy file names.

    Specifications:
    https://msdn.microsoft.com/en-us/library/windows/desktop/bb776902(v=vs.85).aspx#CF_HDROP
    """

    file_bytes = [str.encode(f) for f in file_names]
    file_name_buffer = str.encode('\0').join(file_bytes) + str.encode('\0\0')
    fmt = "lllll%ss" % len(file_name_buffer)
    dropfiles = struct.pack(fmt, 20, 0, 0, 0, 0, file_name_buffer)
    return dropfiles


def copy_file_to_clipboard(file_name):

    def enum_windows_callback(handle, handles):
        if MY_WINDOW_NAME in win32gui.GetWindowText(handle):
            handles.append(handle)

    # Get the handle of my window
    handles = []
    win32gui.EnumWindows(enum_windows_callback, handles)
    my_handle = handles[0]

    while True:
        try:
            # Store the foreground
            foreground_handle = win32gui.GetForegroundWindow()
            # Make sure we can change the foreground
            foreground_tid = win32process.GetWindowThreadProcessId(foreground_handle)[0]
            my_tid = win32api.GetCurrentThreadId()

            win32process.AttachThreadInput(foreground_tid, my_tid, win32con.TRUE)
            ctypes.windll.user32.AllowSetForegroundWindow(-1)
            win32gui.SystemParametersInfo(win32con.SPI_SETFOREGROUNDLOCKTIMEOUT, 0, win32con.SPIF_SENDWININICHANGE | win32con.SPIF_UPDATEINIFILE)
            # Set it
            win32gui.SetForegroundWindow(my_handle)
            win32process.AttachThreadInput(foreground_tid, my_tid, win32con.FALSE)
            break
        except:
            print("Did not work, trying again...")

    # Copy to the clipboard
    file_name = file_name.replace('/', '\\')
    win32clipboard.OpenClipboard(0)
    win32clipboard.EmptyClipboard()
    win32clipboard.SetClipboardData(win32con.CF_HDROP, to_cf_hdrop((file_name,)))
    win32clipboard.CloseClipboard()

    # Change the foreground back
    win32gui.SetForegroundWindow(foreground_handle)

class Something:

    def __init__(self):
        self.d = dispatcher.Dispatcher()

        # Map channels to functions
        self.d.map("/read_sample", self.read_sample)

        # Begin listening on the specified port
        self.server = osc_server.ThreadingOSCUDPServer((IP, PORT_RECEIVE), self.d)

    def run(self):
        self.server.serve_forever()

    def read_sample(self, channel, input_sample):
        copy_file_to_clipboard(input_sample)

        # client = udp_client.SimpleUDPClient(IP, 5005)

        # for x in range(10):
            # time.sleep(1)

if __name__ == "__main__":
    # x = 'A:\\User Library\\SampleLibrary\\Airhead\\I will stay away.alc'
    # copy_file_to_clipboard(x)
    Something().run()
    # print("Sleeping...")

    # def get_my_handle():


        # windows = []
        # pid = os.getpid()
        # print("My Pid: ", pid)
        # win32gui.EnumWindows(callback, (pid, windows))
        # print(windows)
        # return windows[0]

    # foreground_handle = win32gui.GetForegroundWindow()
    # my_handle = get_my_handle()
    # asdf
    
    # win32gui.SetForegroundWindow(self.handle)
    # # Do shit
    # win32gui.SetForegroundWindow(my_handle)
    
    # import time
    # time.sleep(5)

    # print("Done")
    # y = 'A:\\User Library\\SampleLibrary\\Airhead\\Spooky vocals and dropping.alc'
    # copy_file_to_clipboard(y)

