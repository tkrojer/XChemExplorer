FROM reskyner/ccp4
#RUN apt-get install -f libxcursor-dev
ADD ./run_tests .
ADD ./compile_test.py .
RUN git clone https://github.com/reskyner/XChemExplorer.git
#CMD ccp4-python ./compile_tests.py
