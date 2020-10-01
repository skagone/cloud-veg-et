# Gold Directory to be replicated/cloned/cp
- This is the Gold Directory for Creating Custom Docker Runs for ET

[MAGIC](#lots-of-magic)
# TL;DR

```
make model_prepare
```

```
python3 run_drb_aoi.py
bash cmd_runner_drb_tile_7_20.sh

```

# What does it really do?

## Lots of Magic

1. Creates a custom docker image with the **latest** software CLONED from skagone/cloud_veg_et
2. Divides the run into smaller geographic 2 degree BY 2 degree chunks - hardcoded for DRB!
3. creates a BASH shell script file with the complicated docker run commands **(TIMES 12)** for DRB.

## Safety Tip

- tmux

> use tmux to maintain your session over a VPN (with timeouts)

## The code
- gridmeister.py is the main worker for this magic
  - its modular so don't be scared!

### See the Makefile for more info

# BUGS and Warnings

 - for logs to actually see the outputs use bash cmd_runner.......sh
 - used to be that the tile must have the string 'tile' for this software to work - not relevant here

