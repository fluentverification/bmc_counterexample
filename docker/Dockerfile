# syntax=docker/dockerfile:1

FROM movesrwth/stormpy:stable
WORKDIR /opt
RUN pip install z3-solver
RUN mkdir repo
WORKDIR /opt/repo
RUN git clone --depth 1 --branch qest2024 https://github.com/fluentverification/bmc_counterexample.git
WORKDIR /opt/repo/bmc_counterexample
