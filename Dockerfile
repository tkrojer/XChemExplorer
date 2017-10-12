FROM reskyner/ccp4
RUN apt-get install libxcursor-dev
ADD ./run_tests
