# Sample Sorter

## Dependencies

### Windows

Make sure you have [Python 3](https://www.python.org/downloads/windows/) installed.
Then install the following python packages via pip.

    pip install tensorflow python-osc pywin32

Fetch the Windows library from

## Running

Launch the python server in the back

    cd SampleSorter/python_src
    python to_osc.py

Then open Ableton Live and drag the Max 4 Live patch into a new track.

## Compiling from source

Compiling the C++ library requires the following libraries:

- boost
- tinyxml2
- fftw3

On Linux or OS X get the required dependencies from your standard package manager.
On Windows get the required dependencies from ```vcpkg```.

    cd build
    cmake -G "Visual Studio 15 2017 Win64" .. --DCMAKE_TOOLDCHAIN_FILE=A:\vcpkg\scripts\buildsystems\vcpkg.cmake
    cmake --build . --config Releasse
