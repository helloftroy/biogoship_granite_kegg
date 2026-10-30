"""Every GPU joins one NCCL process group, sums a tensor, and times a 1 GB all-reduce.

    torchrun --standalone --nproc_per_node=4 cluster/ptolemy/ddp_check.py
"""

import os
import time

import torch
import torch.distributed as dist


def main() -> None:
    dist.init_process_group("nccl")
    rank, world = dist.get_rank(), dist.get_world_size()
    torch.cuda.set_device(int(os.environ["LOCAL_RANK"]))
    print(f"rank {rank}/{world}: {torch.cuda.get_device_name()}", flush=True)

    x = torch.full((1,), float(rank), device="cuda")
    dist.all_reduce(x)
    expected = world * (world - 1) / 2
    big = torch.ones(256 * 1024 * 1024, device="cuda")  # 1 GB of float32
    dist.all_reduce(big)  # warm-up
    torch.cuda.synchronize()
    start = time.perf_counter()
    for _ in range(5):
        dist.all_reduce(big)
    torch.cuda.synchronize()
    seconds = (time.perf_counter() - start) / 5
    if rank == 0:
        bus = 2 * (world - 1) / world * big.numel() * 4 / seconds / 1e9  # NCCL "bus bandwidth"
        print(f"sum check: {x.item():g} (expected {expected:g}) {'OK' if x.item() == expected else 'WRONG'}")
        print(f"1 GB all-reduce: {seconds * 1000:.0f} ms, bus bandwidth {bus:.0f} GB/s (NVLink A100: ~200+; PCIe: ~10-20)")
    dist.destroy_process_group()


if __name__ == "__main__":
    main()
