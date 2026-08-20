# Portability design

Maintained code follows these rules:

- no `/home/<user>/...`, `~/data/...`, or `/mnt/data/...` paths;
- no hard-coded VAMPIRE installation path;
- all input/output paths are CLI arguments or repository-relative defaults;
- VAMPIRE is resolved with `--vampire`, `VAMPIRE_BIN`, then PATH;
- generated data are written to caller-selected locations and ignored by git;
- one-off convergence/debugging scripts remain under `archive/` and are not advertised as production code.

The only project-specific defaults that remain are *scientific defaults* (for example the thesis model files,
CMC step counts, and the calibrated torque-density volume). They are exposed as command-line options where
changing them is physically meaningful.
