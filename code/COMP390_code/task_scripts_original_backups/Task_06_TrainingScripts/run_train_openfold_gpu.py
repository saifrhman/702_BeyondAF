import os
import sys
import runpy

import pytorch_lightning as pl

# --- 强制 Trainer 默认用 GPU ---
_OrigTrainer = pl.Trainer

def _Trainer(*args, **kwargs):
    # Lightning 2.x: accelerator/devices
    kwargs.setdefault("accelerator", "gpu")
    kwargs.setdefault("devices", 1)

    # 如果你想让它跟着 CUDA_VISIBLE_DEVICES 自动走，也可以不写 devices
    # kwargs.setdefault("devices", "auto")

    return _OrigTrainer(*args, **kwargs)

pl.Trainer = _Trainer

# 可选：避免 triton 缓存在 NFS（性能/退出卡住）
# os.environ.setdefault("TRITON_CACHE_DIR", "/users/sgmwu14/fastscratch/.triton_cache")

# --- 直接以 __main__ 方式执行原脚本，不改一行源码 ---
runpy.run_path("/users/sgmwu14/openfold/train_openfold.py", run_name="__main__")