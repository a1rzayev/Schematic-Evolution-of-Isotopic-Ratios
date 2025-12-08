# Schematic-Evolution-of-Isotopic-Ratios

## What it is
- A tiny desktop window (Tkinter) that shows two plots: Hf and Nd isotopic evolution.
- You set “Crust extraction time (Ga)” and click “Generate plots”. Switch Hf/Nd with the radio buttons. “Save PNGs” writes both images.

## Quick start (macOS, system Python)
```
python3 Nd-Sr_Anticorrelation.py
```

## Quick start (Windows, python.org install)
```
cd C:\path_to\Schematic-Evolution-of-Isotopic-Ratios
python -m venv .venv
.\.venv\Scripts\activate
pip install numpy matplotlib
python Nd-Sr_Anticorrelation.py
```

## If using the project venv on macOS (needs network once)
```
cd /path_to/Schematic-Evolution-of-Isotopic-Ratios
./.venv/bin/python -m pip install numpy matplotlib
MPLCONFIGDIR=$PWD/.mplconfig XDG_CACHE_HOME=$PWD/.cache ./.venv/bin/python Nd-Sr_Anticorrelation.py
```

## Tips
- If the window is empty/gray, check the terminal for errors; Tk must be available.
- The generated files are `evolution_Hf.png` and `evolution_Nd.png`.