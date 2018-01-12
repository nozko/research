#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
from datetime import datetime
import argparse
import sys

import road_simulation


class control:

	# control from snowfall sensor
	def judge_1(self, pre):
		if(pre > 0):		heater = 1
		else:				heater = 0
		return heater


	# control from snow accumulation
	def judge_2(self, snow):
		if(snow > 0.05):	heater = 1
		else:				heater = 0
		return heater


	# always off
	def off(self):
		heater = 0
		return heater



class Qlearning:
	# parameters
	alpha = 0.1
	gamma = 0.9

	def __init__(self):
		self.action = ['on', 'off']

		# rewards
		self.r_on = -1
		self.r_no = 100


	def initializeQ(self, mode):
		# only snow accumulation
		if(mode==0):
			Q = np.zeros(11)

		# snow accumulation and temperature
		elif(mode==1):
			Q = np.zeros((11, 15))

		# invalid mode
		else:
			print('invalid mode error!')
			sys.exit()


	def calcQ(self):
		pass
