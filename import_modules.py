import os, sys, glob
from datetime import datetime
import time

sys.path.append(os.path.join(os.getenv('XChemExplorer_DIR'),'lib'))
sys.path.append(os.path.join(os.getenv('XChemExplorer_DIR'),'web'))
sys.path.append(os.path.join(os.getenv('XChemExplorer_DIR'),'gui'))
sys.path.append(os.path.join(os.getenv('XChemExplorer_DIR'),'gui/maintab_funcs'))

import subprocess
import pickle
import base64
import math
import multiprocessing
import webbrowser
import getpass

from PyQt4 import QtGui, QtCore, QtWebKit
from XChemUtils import parse
from XChemUtils import external_software
from XChemUtils import helpers
import XChemThread
import XChemDB
import XChemPANDDA
import XChemToolTips
import XChemMain
import XChemPlots
import XChemLog
import XChemProcess
import XChemDeposit
import XChemWeb
import dataset_functions as datafunc

import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
import numpy as np
