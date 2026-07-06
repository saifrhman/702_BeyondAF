import os
import runpy
from pathlib import Path

import pytorch_lightning as pl

# Force PyTorch Lightning Trainer to use one GPU by default.
_OrigTrainer = pl.Trainer


def _Trainer(*args, **kwargs):
    kwargs.setdefault("accelerator", "gpu")
    kwargs.setdefault("devices", 1)
    return _OrigTrainer(*args, **kwargs)


pl.Trainer = _Trainer

home = Path.home()
comp702_root = Path(os.environ.get("COMP702_ROOT", home / "COMP702_BeyondAF"))
comp390_root = Path(os.environ.get("COMP390_ROOT", comp702_root / "code" / "COMP390_code"))
openfold_code_dir = Path(os.environ.get("OPENFOLD_CODE_DIR", comp390_root / "openfold"))

triton_cache_dir = os.environ.get(
    "TRITON_CACHE_DIR",
    str(home / "fastscratch" / "triton_cache"),
)
os.environ.setdefault("TRITON_CACHE_DIR", triton_cache_dir)

train_script = openfold_code_dir / "train_openfold.py"

if not train_script.is_file():
    raise FileNotFoundError(f"OpenFold training script not found: {train_script}")

runpy.run_path(str(train_script), run_name="__main__")
