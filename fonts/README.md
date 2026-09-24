# Helvetica figure fonts

All ordinary figure text uses genuine Helvetica. Equations retain their existing
math renderer. Helvetica is a separately supplied font, not a Conda dependency.

Install Helvetica on your system, or copy your Helvetica `.ttf`/`.otf` font files
into `fonts/Helvetica/`. Include regular, bold, oblique/italic, and bold oblique/
italic styles to preserve the figure emphasis. The supported internal font
family names are `Helvetica` and `Helvetica LT Std`.

The figure scripts register files from this directory at startup, so a system-wide
font installation is not required. Restart Python after adding the files, then run:

```bash
conda activate sep_2026
make check-env
make figures
```

Font binaries are ignored by Git and must be supplied by each person regenerating
the figures. The scripts stop with an actionable error if Helvetica is missing;
Nimbus Sans, Arial, and DejaVu Sans are not used as substitutes.
