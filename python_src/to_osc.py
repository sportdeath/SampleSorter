import tensorflow as tf
import numpy as np

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

from sample_reader import SampleReader
from octave_classifier import OctaveClassifier

IP = "127.0.0.1"
PORT_SEND = 5005
PORT_RECEIVE = 5006
MY_WINDOW_NAME = "Command Prompt"
USER_LIBRARY = "A:/User Library/"
SAMPLE_LIBRARY = "A:/User Library/SampleLibrary/"
FORCE_REPROCESS = False
TENSORFLOW_GRAPH = "graph"
TUNING_EPSILON = 15
MOST_SEMITONES = 3
DECISION_THRESHOLD = 0.7
PAGE_NUM = 6

def filter_samples(octaves, tempos, tunings, paths, audio_paths, input_tempo):
    # Determine if tuning is ~semitone
    delta_tunings = 1200. * np.log2(input_tempo/tempos)
    delta_tunings = np.mod(delta_tunings + 600., 1200.) - 600.
    is_in_tune = (50. - np.abs(50. - np.mod(delta_tunings, 100.))) < TUNING_EPSILON
    is_in_tune = np.logical_and(is_in_tune, (np.abs(delta_tunings) < 100 * MOST_SEMITONES))

    octaves = octaves[is_in_tune]
    tempos = tempos[is_in_tune]
    tunings = tunings[is_in_tune]
    delta_tunings = delta_tunings[is_in_tune]
    paths = [path for (path, b) in zip(paths, is_in_tune) if b]
    audio_paths = [path for (path, b) in zip(audio_paths, is_in_tune) if b]

    # Roll the octave
    discrete_tunings = np.round(delta_tunings/100).astype(int)

    rows, cols = np.ogrid[:octaves.shape[0], :octaves.shape[1]]

    discrete_tunings = np.mod(discrete_tunings, octaves.shape[1])
    cols = cols - discrete_tunings[:,np.newaxis]

    octaves = octaves[rows, cols]

    tunings = tunings + delta_tunings

    return octaves, tempos, tunings, paths, audio_paths

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
        print("Loading sample library...")
        self.octaves, self.tempos, self.tunings, self.paths, self.audio_paths = SampleReader.read_samples(
                USER_LIBRARY,
                SAMPLE_LIBRARY,
                FORCE_REPROCESS,
                shuffle=False)
        print(self.octaves.shape[0], "samples loaded.")

        # Load the tensorflow graph from the saved state
        # self.octave_ph, self.training_ph, self.decision = OctaveClassifier().construct()
        saver = tf.train.import_meta_graph(TENSORFLOW_GRAPH + ".meta")
        self.sess = tf.Session()
        saver.restore(self.sess, TENSORFLOW_GRAPH)

        # Get the inputs and outputs of the graph
        graph = tf.get_default_graph()
        self.octaves_ph = graph.get_tensor_by_name("octave_classifier/octave:0")
        self.training_ph = graph.get_tensor_by_name("octave_classifier/training:0")
        self.decision = graph.get_tensor_by_name("octave_classifier/decision:0")

        # Set up an OSC dispatcher
        self.d = dispatcher.Dispatcher()

        # map channels to functions
        self.d.map("/read_sample", self.read_sample)
        self.d.map("/copy_sample", self.copy_sample)
        self.d.map("/load_page", self.load_page)

        # Begin listening on the specified port
        self.server = osc_server.ThreadingOSCUDPServer((IP, PORT_RECEIVE), self.d)

        # Create a client
        self.client = udp_client.SimpleUDPClient(IP, PORT_SEND)

    def run(self):
        self.server.serve_forever()

    def copy_sample(self, channel, sample_index):
        copy_file_to_clipboard(self.samples[sample_index][2])

    def read_sample(self, channel, input_sample):
        octave, tempo, tuning, audio_path = SampleReader.read_sample(input_sample, USER_LIBRARY, FORCE_REPROCESS)

        octave = octave/np.linalg.norm(octave)
        tempo = tempo

        self.client.send_message("/octave", [x/max(octave) for x in octave])
        self.client.send_message("/tempo", 60. * tempo)
        self.client.send_message("/tuning", tuning)

        octaves, tempos, tunings, paths, audio_paths = filter_samples(self.octaves, self.tempos, self.tunings, self.paths, self.audio_paths, tempo)

        # Combine this sample with the other samples
        input_octaves = octave + octaves
        input_octaves = input_octaves/np.linalg.norm(input_octaves, axis=1)[:, None]

        # Feed the samples into the graph and run it
        feed_dict = {}
        feed_dict[self.octaves_ph] = input_octaves
        feed_dict[self.training_ph] = False
        decision_ = self.sess.run(tf.sigmoid(self.decision), feed_dict=feed_dict)
        decision_ = list(decision_[:,0])

        # Combine the results with paths
        zipped = zip(decision_, tunings, paths, audio_paths)
        self.samples = sorted(zipped, key=lambda x: x[0])[::-1]

        self.client.send_message("/playlist", "clear")

        self.load_page("", 0)

        copy_file_to_clipboard(input_sample)

    def load_page(self, channel, page_index):
        for i in range(PAGE_NUM):
            index = page_index * PAGE_NUM + i
            decision, tuning, path, audio_path = self.samples[index]
            if decision > DECISION_THRESHOLD:
                self.client.send_message("/play" + str(i), "open \"" + audio_path + "\"")
                self.client.send_message("/play" + str(i), "speed " + str(2. ** (tuning/1200.)))
                self.client.send_message("/name" + str(i), path.replace('\\','/').replace(SAMPLE_LIBRARY, '')[:-4] + "@" + str(tuning))

if __name__ == "__main__":
    Something().run()
