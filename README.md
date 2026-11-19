```bash
echo "deb [trusted=yes] https://raw.githubusercontent.com/christian-rauch/exotica/jammy-humble-amd64/ ./" | sudo tee /etc/apt/sources.list.d/christian-rauch_exotica-jammy-humble-amd64.list
echo "yaml https://github.com/christian-rauch/exotica/raw/jammy-humble-amd64/local.yaml humble" | sudo tee /etc/ros/rosdep/sources.list.d/1-christian-rauch_exotica-jammy-humble-amd64.list
```
