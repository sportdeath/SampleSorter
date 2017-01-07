#ifndef MP4_FILE_H
#define MP4_FILE_H

class MP4File {
  public:
    static std::vector<std::vector<double> >
      read(std::string filePath,
           double startSeconds,
           double endSeconds);
};

#endif
