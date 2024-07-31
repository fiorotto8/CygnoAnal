# Directionality and Other features

Analyzer class and methods useful for directionality and general analysis on CYGNO reconstructed tracks

## To compile

```g++ Analyzer.cxx Base_script.cxx -o nameprog `root-config --libs --cflags` -lSpectrum```

Optimized:
```g++ Analyzer.cxx Base_script.cxx -O3 -o nameprog `root-config --libs --cflags` -lSpectrum```

## To run

```./nameprog path_to_rootfile [output_directory]```

## Python version

- a minimal python implementation of the Analyzer class is in `Analyzer.py`
- example for ALpha diffusion analysis is in `DiffusionAlpha.py` and a script to aggregate and plot in `plotDiffusion.py`