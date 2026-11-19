```bash
echo "deb [trusted=yes] https://raw.githubusercontent.com/christian-rauch/exotica/noble-jazzy-amd64/ ./" | sudo tee /etc/apt/sources.list.d/christian-rauch_exotica-noble-jazzy-amd64.list
echo "yaml https://github.com/christian-rauch/exotica/raw/noble-jazzy-amd64/local.yaml jazzy" | sudo tee /etc/ros/rosdep/sources.list.d/1-christian-rauch_exotica-noble-jazzy-amd64.list
```
