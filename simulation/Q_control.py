#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
from datetime import datetime
import argparse

import road_simulation


class control:

	# control from snowfall sensor
	def judge_1(self, pre):
		if(pre > 0):		heater = 1
		else:				heater = 0
		return heater


	# control from snow accumulation
	def judge_2(self, snow):
		if(snow > 0.03):	heater = 1
		else:				heater = 0
		return heater


	# always off
	def off(self):
		heater = 0
		return heater



class Qlearning:
	alpha = 0.1
	gamma = 0.9

	def __init__(self):
		self.action = ['on', 'off']

		# rewards
		self.r_on = -1
		self.r_no = 100


	def calcQ(self):
		pass
