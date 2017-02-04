import ctypes
import math

lib = ctypes.CDLL("../build/libSampleSorter.so")

class Sample:
    def __init__(self, fileName, userLibrary):
        lib.NewAbletonSampleFile.restype = ctypes.c_void_p
        self.s = lib.NewAbletonSampleFile(str(fileName).encode('ascii'), str(userLibrary).encode('ascii'))

    def getFileName(self):
        lib.getFileName.argtypes = [ctypes.c_void_p]
        lib.getFileName.restype = ctypes.c_char_p
        return lib.getFileName(self.s)
        
    def process(self):
        lib.process.argtypes = [ctypes.c_void_p]
        lib.process.restype = ctypes.c_bool
        return lib.process(self.s)

    def getTuning(self):
        lib.getTuningCents.argtypes = [ctypes.c_void_p]
        lib.getTuningCents.restype = ctypes.c_long
        return lib.getTuningCents(self.s)

    def getTheOne(self):
        lib.getTheOneWithTuning.argtypes = [ctypes.c_void_p]
        lib.getTheOneWithTuning.restype = ctypes.c_double
        return lib.getTheOneWithTuning(self.s)

    def getTempo(self):
        lib.getBeatWithTuning.argtypes = [ctypes.c_void_p]
        lib.getBeatWithTuning.restype = ctypes.c_double
        return lib.getBeatWithTuning(self.s)

    def getChords(self):
        # get the number
        lib.getNumChords.argtypes = [ctypes.c_void_p]
        lib.getNumChords.restype = ctypes.c_size_t
        numChords = lib.getNumChords(self.s)

        # allocate enough space
        lib.getChords.argtypes = [ctypes.c_void_p]
        lib.getChords.restype = ctypes.POINTER(ctypes.POINTER(ctypes.c_double * 12) * numChords)
        c_chords = lib.getChords(self.s)
        chords = [[value for value in chord.contents] for chord in c_chords.contents]

        # delete allocated space
        lib.deleteChords.argtypes = [ctypes.c_void_p, ctypes.c_size_t]
        lib.deleteChords(c_chords, numChords)

        # return
        return chords

    def writeToFile(self):
        lib.writeToFile.argtypes = [ctypes.c_void_p]
        lib.writeToFile(self.s)

    def delete(self):
        lib.deleteAbletonSampleFile.argtypes = [ctypes.c_void_p]
        lib.deleteAbletonSampleFile(self.s)


    def writeToMIDI(self, midiFileName):
        import midi
        chords = self.getChords()
        # find the maximum amplitude
        maximumAmp = 0
        for c in range(len(chords)):
            for b in range(12):
                maximumAmp = max(maximumAmp, math.sqrt(chords[c][b]))

        MIDI_velocity_ratio = 127./maximumAmp
        ticks_per_quarter_note = 220

        # initialize MIDI
        pattern = midi.Pattern()
        track = midi.Track()
        pattern.append(track)

        # add Tempo
        tempo = midi.SetTempoEvent(tick = 0, bpm = self.getTempo() * 60)
        track.append(tempo)

        for c in range(len(chords)):
            # Add note ons
            for b in range(12):
                vel = MIDI_velocity_ratio*math.sqrt(chords[c][b])
                on = midi.NoteOnEvent(tick = 0, velocity = int(vel), pitch=midi.A_4+b)
                track.append(on)

            # Add note offs
            off = midi.NoteOffEvent(tick = ticks_per_quarter_note, pitch=midi.A_4)
            track.append(off)
            for b in range(1,12):
                off = midi.NoteOffEvent(tick = 0, pitch=midi.A_4+b)
                track.append(off)

        # Add the end of the file
        eot = midi.EndOfTrackEvent(tick=1)
        track.append(eot)

        # write to file
        midi.write_midifile(midiFileName, pattern)
