################## Base Image ##########
ARG PYTHON_VERSION="3.10"
FROM --platform=linux/amd64 python:${PYTHON_VERSION}

################## ARGUMENTS/Environments ##########
ARG BUILD_DATE
ARG BUILD_VERSION
ARG LICENSE="Apache-2.0"
ARG iAnnotateSV_VERSION
ARG VCS_REF

################## METADATA ########################
LABEL org.opencontainers.image.vendor="MSKCC"
LABEL org.opencontainers.image.authors="Ronak Shah (shahr2@mskcc.org)"

LABEL org.opencontainers.image.created=${BUILD_DATE} \
    org.opencontainers.image.version=${BUILD_VERSION} \
    org.opencontainers.image.licenses=${LICENSE} \
    org.opencontainers.image.version.pvs=${iAnnotateSV_VERSION} \
    org.opencontainers.image.source.pv="https://pypi.org/project/iAnnotateSV" \
    org.opencontainers.image.vcs-url="https://github.com/rhshah/iAnnotateSV.git" \
    org.opencontainers.image.vcs-ref=${VCS_REF}

LABEL org.opencontainers.image.description="This container uses python 3.10 as the base image to build \
    iAnnotateSV version ${iAnnotateSV_VERSION}"

################## INSTALL ##########################

WORKDIR /app
ADD . /app

# install iAnnotateSV 

RUN apt-get update && apt-get install --no-install-recommends -y gcc g++ zlib1g-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && make deps-install \ 
    && poetry build \
    && pip install dist/iAnnotateSV-*-py3-none-any.whl
