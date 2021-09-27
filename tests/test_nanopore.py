#!/usr/bin/env python
import unittest
import logging
import filecmp
import mock
import os

from taca.nanopore.nanopore import *
from taca.utils import config as conf