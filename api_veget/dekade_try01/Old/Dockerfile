FROM python

ENV VERS=RUNN

RUN apt-get update && \
	apt-get install -y vim

RUN \
	pip install --no-cache matplotlib && \
	pip install --no-cache matplotlib_venn && \
	pip install --no-cache pandas && \
	pip install --no-cache seaborn && \
	pip install --no-cache plotly && \
	pip install --no-cache plotly.express && \
	pip install --no-cache rasterio && \
	pip install --no-cache awscli && \
	pip install --no-cache s3fs && \
	pip install --no-cache ffmpeg && \
	pip install --no-cache pika && \
	pip install --no-cache dash && \
	pip install --no-cache boto3 && \
	pip install --no-cache fiona


ENV TONYVERS=6.9
RUN mkdir -p /home/veget \
	&& cd /home/veget\
	&& git clone http://github.com/skagone/cloud-veg-et \
	&& (cd /home/veget/cloud-veg-et/vegetLib; make build-root) \
	&& mkdir -p /home/veget/cloud-veg-et/api_veget/log

COPY running_config /home/veget/cloud-veg-et/api_veget/running_config

RUN apt-get install ca-certificates && mkdir -p /etc/pki/tls/certs && \
	cp /etc/ssl/certs/ca-certificates.crt /etc/pki/tls/certs/ca-bundle.crt

WORKDIR /home/veget/cloud-veg-et/api_veget
