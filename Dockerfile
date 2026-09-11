# A pinned environment for reproducing the shipped example records.
#
# The package itself needs only the ranges in pyproject.toml. This image
# instead pins the exact versions the shipped records were computed with (the
# release stack of REPRODUCIBILITY.md), so that a difference between a rerun
# and the shipped record is a difference in the work and not in the stack. The two are separate claims and the file carrying each
# is separate as well: pyproject.toml states what the software supports,
# requirements-lock.txt states the pinned stack the numbers reproduce on.
#
# Build and check:
#   docker build -t amr-clonalshare:1.0.0 .
#   docker run --rm amr-clonalshare:1.0.0 --version
#   docker run --rm amr-clonalshare:1.0.0 \
#       --config examples/ssuis/config.yaml --results-dir /tmp/ssuis
# Your own data:
#   docker run --rm -v "$PWD":/data amr-clonalshare:1.0.0 \
#       --config /data/config.yaml --results-dir /data/out
FROM python:3.11.15-slim-bookworm

LABEL org.opencontainers.image.title="amr-clonalshare"
LABEL org.opencontainers.image.version="1.0.0"
LABEL org.opencontainers.image.licenses="MIT"
LABEL org.opencontainers.image.description="The share of an antimicrobial resistance phenotype that the lineages of a bacterial collection carry; pinned stack for the article"

ENV PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1 \
    OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

WORKDIR /opt/amr-clonalshare
COPY requirements-lock.txt ./
RUN python -m pip install --no-cache-dir --upgrade pip \
 && python -m pip install --no-cache-dir -r requirements-lock.txt

COPY pyproject.toml README.md LICENSE Licence.txt MANIFEST.in ./
COPY src ./src
COPY examples ./examples
RUN python -m pip install --no-cache-dir --no-deps . \
 && useradd --create-home --uid 1000 runner \
 && chown -R runner:runner /opt/amr-clonalshare
USER runner

ENTRYPOINT ["amr-clonalshare"]
CMD ["--help"]
