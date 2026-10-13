"""Load the IBM Granite Geospatial Ocean encoder and turn tiles into embeddings.

The checkpoint is a masked autoencoder; only its encoder is used. It is built as a
TerraTorch `prithvi_eo_tiny` backbone with the Granite configuration (17 channels,
42x42 input, 2x2 patches, 512-wide, 12 layers) and the encoder weights are loaded
from the Hugging Face snapshot. An embedding is the mean of the encoder's output
tokens: one 512-vector per tile.

Requires torch, terratorch and huggingface_hub.
"""

from __future__ import annotations

import importlib
import os
from pathlib import Path

import numpy as np
import torch
from torch import nn

MODEL_NAME = "ibm-granite/granite-geospatial-ocean"
BACKBONE_NAME = "prithvi_eo_tiny"
EMBED_DIM = 512
# Parts of the autoencoder checkpoint that are not the encoder.
DECODER_MARKERS = ("decoder", "mask_token", "pred", "head", "norm_pix_loss")
KEY_PREFIXES = ("module.", "model.", "backbone.", "encoder.")


def build_encoder(model_name: str = MODEL_NAME, local_files_only: bool = False) -> tuple[nn.Module, str]:
    """The Granite encoder with pretrained weights, plus a one-line report of what was loaded."""
    from huggingface_hub import snapshot_download

    try:
        snapshot_dir = Path(snapshot_download(model_name, local_files_only=local_files_only))
    except Exception as exc:
        if local_files_only:
            raise RuntimeError(
                f"{model_name} is not in the local Hugging Face cache (HF_HOME={os.environ.get('HF_HOME', '<unset>')}). "
                "Run cluster/cache_granite.sbatch on an internet-enabled node first."
            ) from exc
        raise
    checkpoint_path = snapshot_dir / "checkpoint.pt"
    if not checkpoint_path.exists():
        raise FileNotFoundError(f"Expected Granite checkpoint not found: {checkpoint_path}")

    encoder = _backbone_registry().build(
        BACKBONE_NAME,
        pretrained=False,
        bands=list(range(17)),
        in_chans=17,
        num_frames=1,
        img_size=(42, 42),
        patch_size=(1, 2, 2),
        embed_dim=EMBED_DIM,
        depth=12,
        num_heads=8,
        mlp_ratio=4,
    )
    checkpoint = torch.load(checkpoint_path, map_location="cpu")
    loadable, skipped = match_encoder_weights(_state_dict(checkpoint), encoder.state_dict())
    if not loadable:
        raise RuntimeError("No checkpoint tensors matched the encoder; check key prefixes and configuration.")
    result = encoder.load_state_dict(loadable, strict=False)
    not_decoder = [key for key, reason in skipped if reason != "decoder"]
    report = (
        f"Loaded Granite encoder from {checkpoint_path}: {len(loadable):,} tensors loaded, "
        f"{len(skipped) - len(not_decoder):,} decoder tensors skipped, {len(not_decoder):,} other tensors unmatched, "
        f"{len(result.missing_keys):,} encoder tensors left at initial values"
    )
    if not_decoder or result.missing_keys:
        report += "\n  unmatched checkpoint keys: " + ", ".join(not_decoder[:10])
        report += "\n  encoder keys without weights: " + ", ".join(list(result.missing_keys)[:10])
    return encoder, report


def _backbone_registry():
    try:
        importlib.import_module("terratorch.models.backbones")
    except ImportError:
        pass
    errors = []
    for module_name in (
        "terratorch.registry",
        "terratorch.models.registry",
        "terratorch.models.backbones",
        "terratorch.models.backbones.registry",
    ):
        try:
            registry = getattr(importlib.import_module(module_name), "BACKBONE_REGISTRY", None)
        except ImportError as exc:
            errors.append(f"{module_name}: {exc}")
            continue
        if registry is not None:
            return registry
    raise RuntimeError("Could not import TerraTorch's BACKBONE_REGISTRY (pip install terratorch).\n" + "\n".join(errors))


def _state_dict(checkpoint) -> dict[str, torch.Tensor]:
    if isinstance(checkpoint, dict):
        for key in ("state_dict", "model", "model_state_dict", "checkpoint"):
            if isinstance(checkpoint.get(key), dict):
                return checkpoint[key]
        if all(isinstance(value, torch.Tensor) for value in checkpoint.values()):
            return checkpoint
    raise ValueError("Could not find a tensor state_dict inside checkpoint.pt")


def match_encoder_weights(
    raw_state: dict[str, torch.Tensor], model_state: dict[str, torch.Tensor]
) -> tuple[dict[str, torch.Tensor], list[tuple[str, str]]]:
    """Map checkpoint tensors onto encoder keys, trying the key with wrapper prefixes stripped.

    Returns (weights to load, [(skipped checkpoint key, reason)]).
    """
    loadable, skipped = {}, []
    for raw_key, value in raw_state.items():
        if not isinstance(value, torch.Tensor):
            skipped.append((raw_key, "not a tensor"))
            continue
        if any(marker in raw_key for marker in DECODER_MARKERS):
            skipped.append((raw_key, "decoder"))
            continue
        variants, queue = [], [raw_key]
        while queue:
            key = queue.pop(0)
            if key not in variants:
                variants.append(key)
            queue += [key[len(prefix) :] for prefix in KEY_PREFIXES if key.startswith(prefix)]
        match = next((k for k in variants if k in model_state and model_state[k].shape == value.shape), None)
        if match is None:
            skipped.append((raw_key, "no encoder key of the same name and shape"))
        else:
            loadable[match] = value
    return loadable, skipped


def encode(encoder: nn.Module, pixel_values: torch.Tensor) -> torch.Tensor:
    """(batch, 17, 42, 42) tiles -> (batch, 512) embeddings: the mean over the encoder's output tokens."""
    frames = pixel_values.unsqueeze(2)  # add the single-frame time axis the backbone expects
    try:
        outputs = encoder(frames)
    except TypeError:
        try:
            outputs = encoder(x=frames)
        except TypeError:
            outputs = encoder(pixel_values=frames)
    return pool(outputs)


def pool(outputs) -> torch.Tensor:
    if isinstance(outputs, dict):
        for key in ("last_hidden_state", "features", "x", "out"):
            if key in outputs:
                outputs = outputs[key]
                break
        else:
            outputs = next(value for value in outputs.values() if isinstance(value, torch.Tensor))
    if getattr(outputs, "pooler_output", None) is not None:
        return outputs.pooler_output
    if getattr(outputs, "last_hidden_state", None) is not None:
        x = outputs.last_hidden_state
    elif isinstance(outputs, (tuple, list)):
        x = outputs[-1]  # the final layer
    else:
        x = outputs
    if x.ndim == 2:
        return x
    if x.ndim == 3:  # (batch, tokens, dim)
        return x.mean(dim=1)
    if x.ndim in (4, 5):  # (batch, dim, [frames,] height, width)
        return x.mean(dim=tuple(range(2, x.ndim)))
    raise RuntimeError(f"Do not know how to pool encoder output of shape {tuple(x.shape)}")


def load_tiles(paths: list[Path], nan_fill: float = 0.0) -> torch.Tensor:
    """Stack tile .npz files into a batch. NaN (missing SST) becomes nan_fill, i.e. the channel mean."""
    arrays = []
    for path in paths:
        with np.load(path, allow_pickle=False) as tile:
            data = tile["data"].astype(np.float32)
        arrays.append(np.nan_to_num(data, nan=nan_fill, posinf=nan_fill, neginf=nan_fill))
    return torch.from_numpy(np.stack(arrays, axis=0))
