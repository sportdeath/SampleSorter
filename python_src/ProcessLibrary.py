import os

from sample import Sample

library = "/Users/tfh/Dropbox (MIT)/UserLibrary/SampleLibrary/Chops/"

samples = []

for root, directories, filenames in os.walk(library):
    for filename in filenames:
        if filename[-4:] == ".alc":
            filePath = os.path.join(root, filename)
            print "Processing:", filename[:-4]
            s = Sample(filePath)
            sample.append(s)
            print "Tuning:", s.getTuning()
            print "Tempo:", s.getTempo() * 60.
            print

ratios = [1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 3/2.]

for s1 in samples:
    tempo1 = s1.getTempo()
    for s2 in samples:
        tempo2 = s2.getTempo()

        larger = max(tempo1, tempo2)
        smaller = min(tempo1, tempo2)

        ratio = large/small

        for r in ratios:
            # if pitch change < an octave its cool
             



        # round to 2, 3, 4, 6, 8, 12, 16, 24, 32 or 3/2, 4/3, 9/2, 9/4

        # get larger of the two
        # what is whole number ratio large/small
        # how far away is it from that ratio
        # also take into account ratios like 2/3
        
        
