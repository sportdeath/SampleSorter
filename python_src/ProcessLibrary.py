import os, math, signal

from sample import Sample

#import tensorflow as tf

#def getRandomChords(samples):

#def getRandomCHords(pairs):

# to generate a batch:
# for i in batch size:
# choose randomly between regular (1) and combined (-1)
# if randomly choose the intersection point
# randomly rotate the chords
# add and average energy...
# if combined randomly choose start points (its recurrent so train until one of them ends)
# train the output to be 1 consistent
# choose start point randomly and train to the end.

# listen to midi to see if euclidean distance is better than
# simply adding up...

#lstm = tf.nn.rnn_cell.BasicLSTMCELL(lstm_size)
#fw_cell = tf.nn.rnn_cell.MultiRNNCell([lstm] * number_of_layers)
#bw_cell = fw_cell

#initial_state = state = stacked_lstm.zero_state(None, tf.float32)



#x = Placeholder for list of tensors of size [None, 12]

#tf.nn.bidirectional_rnn(fw_cell, bw_cell, x)

def processFiles(library):
    processedFiles = []
    unProcessedFiles = []
    for root, directories, filenames in os.walk(library):
        for filename in filenames:
              if filename[-4:] == ".alc":
                  filePath = os.path.join(root, filename)
                  s = Sample(filePath)
                  if s.process():
                      print "Processed:", s.getFileName()
                      print "Tuning:", s.getTuning()
                      print "Tempo:", s.getTempo() * 60.
                      print
                      s.writeToFile()
                      s.delete()
                      processedFiles.append(filePath)
                  else:
                      print "Could not process ", filename
                      print
    return processedFiles

def findPairs(files, semitonesLowerBound, semitonesUpperBound, tuningBound):
    ratios = [1, 2, 3, 4, 6, 8, 12, 16, 24, 32]
    inverseRatios = [1./r for r in ratios]
    ratios += inverseRatios

    print ratios

    pairs = {}

    for i in range(len(files)):
        s1 = Sample(files[i])
        s1.process()
        tempo1 = s1.getTempo()

        pairs[s1.getFileName()] = []

        for j in range(i+1, len(files)):
            s2 = Sample(files[j])
            s2.process()
            tempo2 = s2.getTempo()

            ratio = tempo1/tempo2

            for r in ratios:
                numCents = 1200. * math.log(ratio/r, 2.)
                distanceFromSemitone = (numCents % 100.)
                if distanceFromSemitone > 50:
                    distanceFromSemitone = 100 - distanceFromSemitone

                if semitonesLowerBound * 100. < numCents < semitonesUpperBound*100. \
                    and distanceFromSemitone < tuningBound:
                    pairs[s1.getFileName()].append((s2.getFileName(), r, numCents))
            s2.delete()
        s1.delete()

    return pairs

def makeBatch(batchSize, files, pairs):
  batch = []
  for i in range(batchSize):
    # choose randomly between file of pair
    # if file:
    # choose a random file 
    # get chords as vector
    # if pair
    # choose a random pair
    # get chords for both
    # choose a random intersection point
    # add them up and divide by 2
    # for both:
    # randomly rotate
    pass
  return batch


def main():
    library = "/Users/tfh/Dropbox (MIT)/UserLibrary/SampleLibrary/Chops/"
    #library = "/Users/tfh/Dropbox (MIT)/SampleSorter/TestFiles/"

    processedFiles = processFiles(library)

    pairs = findPairs(processedFiles, -6, 3, 5)

    print len(pairs),"pairs!"

    maxPairs = 0
    for f in pairs:
      pairsF = len(pairs[f])
      if pairsF > maxPairs:
        maxPairs = pairsF
        maxF = f

    for fileName, ratio, numCents in pairs[maxF]:
      sample = Sample(fileName)
      print sample.getFileName()
      print ratio
      print numCents
      sample.delete()

    samplF = Sample(maxF)
    print samplF.getFileName(), "has", maxPairs,"pairs!"
    samplF.delete()

    


# TODO
# GUI
# find cliques
# drag and drop
# drop files in
# drag files out with appropriate relative tuning
# mp3 support
# warped files support
# straight audio file support
# support for the couple samples that aren't working



if __name__ == "__main__":
  main()
