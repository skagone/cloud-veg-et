Dir=dekade_2
Image=tbutzer/${Dir}

build:
	docker build -t ${Image} .

run:
	docker run -it ${Image} bash

exec:
	docker exec -it ${Image} bash

git:
	(cd /opt/etops ; make)

pull:
	(cd ~/opt/cloud-veg-et; git pull)

kill_all_containers:
	echo WARNING not killing just echo
	echo docker kill '$$(docker ps -q)'

cmd:
	python3 run_drb_aoi.py

config:
	cp ~/dropbox/myconfig/* ./running_config/

model_prepare: config build cmd
	@echo ${Dir}  ${Image}
	@echo TMUX is REALLY Nice -- even when the VPN is Not!
	@echo bash cmd_runner_run_dekade_2.sh
	@echo bash cmd_runner_run_dekade_2.sh



dela_logs:
	docker run -i ${Image} python3 api_veget.py -c running_config -s AOI/d40N-78E.shp tile40N-80E_chip40N-78E


single_tile:
	docker run -i ${Image} python3 api_veget.py -c running_config -s AOI/tile40N-80E_chip32N-72E.shp tile40N-80E_chip32N-72E 


VOLS=-v /opt/etops/grid_veget_c/AOI:/home/veget/cloud-veg-et/veget/api_veget/AOI
mount_volume_test:
	docker run -i ${VOLS} ${Image} python3 api_veget.py -c running_config -s AOI/tile40N-80E_chip32N-72E.shp tile40N-80E_chip32N-72E
