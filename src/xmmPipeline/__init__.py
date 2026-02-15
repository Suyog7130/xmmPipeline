
import os
import sys

# Add directory to path so that Python can find it.
# This also resolves the "module not found" error
# during import from outside the directory (e.g. in Sphinx-autodoc)
# and the error in absolute imports among scripts within directory.

CURR_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(CURR_DIR)

#from . import *
from epicObj import *
from epicPipeline import *
from epicPileup import *
from epicSpectra import *

from findOverlap import *

# TODO: upgrade this to plotutils code, that would be
# made available via `pip install putils`
from plotAnal import *

from convert import *

#__all__ = ['xmmPipeline']
