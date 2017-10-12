FROM reskyner/ccp4
RUN apt-get libxcursor-dev
ADD ./run_tests
