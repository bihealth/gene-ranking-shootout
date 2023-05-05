FROM python:3.8

RUN git clone https://github.com/holtgrewe/CADA.git /app && \
    cd /app && \
    git checkout shootout-fixes && \
    pip install -e .

ENTRYPOINT ["python", "/usr/local/bin/CADA"]
CMD ["--help"]
