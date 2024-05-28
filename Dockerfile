FROM glucksistemi/pymol-ubuntu:2204.2023.11.27
USER root
ADD ./reqirements.txt /app/requirements.txt
RUN pip install -r /app/requirements.txt
COPY pss_worker_framework-0.3.1.tar.gz /pss_worker_framework-0.3.1.tar.gz
RUN pip install /pss_worker_framework-0.3.1.tar.gz
ADD . /app
WORKDIR /app
CMD python runner.py