FROM python:3.11-slim-bookworm

ENV DEBIAN_FRONTEND=noninteractive
WORKDIR /opt/busco2tree

# 1) Paquetes mínimos del sistema (solo lo indispensable para descargar e instalar)
RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates curl bzip2 \
    && rm -rf /var/lib/apt/lists/*

# 2) Instalar micromamba (gestor conda liviano)
#    Lo dejamos en /usr/local/bin para que esté en PATH
RUN curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest \
  | tar -xvj -C /usr/local/bin --strip-components=1 bin/micromamba

# 3) Crear un entorno conda con herramientas bioinfo (bioconda) + libs python
#    - mafft, trimal, iqtree son consistentes en bioconda
#    - biopython, pyyaml los instalamos acá para que queden en el mismo entorno
ENV MAMBA_ROOT_PREFIX=/opt/micromamba
RUN micromamba create -y -n b2t -c conda-forge -c bioconda \
    python=3.11 \
    mafft trimal iqtree \
    biopython pyyaml \
  && micromamba clean -a -y

RUN ln -sf /opt/micromamba/envs/b2t/bin/iqtree3 /usr/local/bin/iqtree

# 4) Asegurar que el entorno quede activo por defecto
ENV PATH=/opt/micromamba/envs/b2t/bin:$PATH

# 5) Copiar tu código
COPY BUSCO2Tree.py ./BUSCO2Tree.py
COPY steps ./steps

# 6) Checks rápidos durante el build (si falla acá, te enterás en build y no en runtime)
RUN python -c "import Bio; import yaml" \
 && mafft --version \
 && trimal -h | head -n 2 \
 && iqtree2 -h | head -n 2

ENTRYPOINT ["python", "/opt/busco2tree/BUSCO2Tree.py"]
