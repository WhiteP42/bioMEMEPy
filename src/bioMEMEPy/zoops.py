import logging
import math
from . import tools
from .tools import BasePWM as PWM

logger = logging.getLogger(__name__)

def e_step():
    raise NotImplementedError

def m_step():
    raise NotImplementedError