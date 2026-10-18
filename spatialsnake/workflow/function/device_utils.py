import torch


def resolve_torch_device(device):
    """Resolve auto/cpu/cuda to a single executable torch device."""
    requested = str(device or "auto").strip().lower()
    if requested not in {"auto", "cpu", "cuda"}:
        raise ValueError("device must be one of: auto, cpu, cuda")

    if requested == "cuda" and not torch.cuda.is_available():
        raise RuntimeError(
            "CUDA was requested, but torch.cuda.is_available() is False. "
            "Use device=auto or device=cpu, or install a compatible CUDA PyTorch build."
        )

    if requested == "auto":
        return "cuda" if torch.cuda.is_available() else "cpu"
    return requested


def lightning_accelerator(device):
    return "gpu" if device == "cuda" else "cpu"
