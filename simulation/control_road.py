#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
from datetime import datetime
import argparse

import road_simulation


class control:

	# snowfall sensor
	def judge_1(self, pre):
		if(pre > 0):		heater = 1
		else:				heater = 0
		return heater


	# control from snow accumulation
	def judge_2(self, snow):
		if(snow > 0):		heater = 1
		else:				heater = 0
		return heater
