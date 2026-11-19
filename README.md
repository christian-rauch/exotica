```bash
echo "deb [trusted=yes] https://raw.githubusercontent.com/christian-rauch/exotica/resolute-lyrical-amd64/ ./" | sudo tee /etc/apt/sources.list.d/christian-rauch_exotica-resolute-lyrical-amd64.list
echo "yaml https://github.com/christian-rauch/exotica/raw/resolute-lyrical-amd64/local.yaml lyrical" | sudo tee /etc/ros/rosdep/sources.list.d/1-christian-rauch_exotica-resolute-lyrical-amd64.list
```
