# Replace the default logging configuration with a custom one
from astropy.logger import logging
log = logging.getLogger('SDTmonitor')
log.propagate = False
sh = logging.StreamHandler()
f = logging.Formatter('%(asctime)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
sh.setFormatter(f)
log.addHandler(sh)
log.setLevel(logging.INFO)

MAX_FEEDS = 7
