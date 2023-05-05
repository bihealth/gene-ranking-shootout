# Gene Ranking Shootout

## Building CADA Podman Image

There is no public REST API or docker image for CADA (yet).
Here is how to build the needed CADA Podman image:

```bash
# cd docker/cada
# bash build.sh
...
# podman run -it --rm localhost/cada-for-shootout:latest --help
```
