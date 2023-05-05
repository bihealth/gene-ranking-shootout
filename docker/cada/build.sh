#!/usr/bin/env bash

set -euo pipefail
set -x

podman build -t cada-for-shootout .

echo "RUN: podman run -it --rm localhost/cada-for-shootout:latest"
