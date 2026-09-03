# A pinned environment for reproducing every number in the manuscript.
#
# The package itself needs only the ranges in pyproject.toml. This image
# instead pins the exact versions the reported results were computed with, so
# that a difference between a rerun and the paper is a difference in the work
# and not in the stack. The two are separate claims and the file carrying each
# is separate as well: pyproject.toml states what the software supports,
# requirements-lock.txt states what produced the numbers.
#
# Build and check:
#   docker build -t amr-clonalshare:1.0.0 .
#   docker run --rm amr-clonalshare:1.0.0 --version
#   docker run --rm amr-clonalshare:1.0.0 \
#       --config examples/synthetic/planted.yaml --results-dir /tmp/planted
# Your own data:
#   docker run --rm -v "$PWD":/data amr-clonalshare:1.0.0 \
#       --config /data/config.yaml --results-dir /data/out
FROM python:3.12.13-slim-bookworm

LABEL org.opencontainers.image.title="amr-clonalshare"
LABEL org.opencontainers.image.version="1.0.0"
LABEL org.opencontainers.image.licenses="MIT"
LABEL org.opencontainers.image.description="Lineage-conditioned diagnostics for bacterial resistance patterns; pinned stack for the SoftwareX manuscript"

ENV PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1 \
    OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

WORKDIR /opt/btc
COPY requirements-lock.txt ./
RUN python -m pip install --no-cache-dir --upgrade pip \
 && python -m pip install --no-cache-dir -r requirements-lock.txt

COPY pyproject.toml README.md LICENSE Licence.txt MANIFEST.in ./
COPY src ./src
COPY examples ./examples
RUN python -m pip install --no-cache-dir --no-deps . \
 && python examples/synthetic/make_data.py \
 && useradd --create-home --uid 1000 btc \
 && chown -R btc:btc /opt/btc
USER btc

ENTRYPOINT ["amr-clonalshare"]
CMD ["--help"]
