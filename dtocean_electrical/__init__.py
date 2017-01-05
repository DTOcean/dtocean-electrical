
import logging

from polite.paths import ObjDirectory, UserDataDirectory, DirectoryMap
from polite.configuration import Logger

# Set default logging handler to avoid "No handler found" warnings.
try:  # Python 2.7+
    from logging import NullHandler
except ImportError:
    class NullHandler(logging.Handler):
        def emit(self, record):
            pass

logging.getLogger(__name__).addHandler(NullHandler())

def start_logging(level=None):

    """Start python logger"""

    objdir = ObjDirectory(__name__, "config")
    datadir = UserDataDirectory("dtocean_electrical", "DTOcean", "config")
    dirmap = DirectoryMap(datadir, objdir)

    log = Logger(dirmap)
    log("dtocean_electrical",
        level=level,
        info_message="Begin logging for dtocean_electrical.")

