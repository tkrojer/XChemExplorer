FROM reskyner/ccp4
#RUN apt-get install -f libxcursor-dev
ADD ./run_tests .
ADD ./compile_test.py .
RUN mkdir XChemExplorer
ADD * /XChemExplorer/
#CMD ccp4-python ./compile_tests.py
