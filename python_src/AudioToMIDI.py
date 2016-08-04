from sample import Sample
import sys

def main(argv):
    inputFile = argv[0]
    outputFile = argv[1]

    s = Sample(inputFile)
    print "Tuning:", s.getTuning()
    print "The One:", s.getTheOne()
    print "Tempo:", s.getTempo() * 60.

    s.writeToMIDI(outputFile)
    s.delete()
    print "written to MIDI"

if __name__ == "__main__":
  main(sys.argv[1:])
