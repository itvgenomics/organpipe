FROM julia:1.11.1

RUN apt-get update && apt-get install -y git

RUN git clone https://github.com/ian-small/chloe /opt/chloe

RUN git clone https://github.com/ian-small/chloe_references /opt/chloe_references

WORKDIR /opt/chloe

ENV JULIA_CPU_TARGET="generic"

ENV JULIA_DEPOT_PATH="/opt/chloe/.julia"

RUN julia --project=. -e 'import Pkg; Pkg.instantiate()'

RUN apt-get update && apt-get install -y python3