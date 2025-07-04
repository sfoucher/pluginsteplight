
To compile the image:

```bash
docker build --network=host -t computree-linux .
```

Execute (tested under WSL2)
```bash
docker run -it --rm --name gcloud --network=host -v /tmp/.X11-unix:/tmp/.X11-unix -e DISPLAY=:0 computree-linux:latest ./Computree.sh
```

You may need to execute `xhost +local:docker`
