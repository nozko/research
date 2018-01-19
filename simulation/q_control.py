#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
from datetime import datetime
import argparse
import sys
import numpy as np
import random

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



action = [0, 1]		# off:0, on:1

class Qlearning:

	def __init__(self, MODE, alpha, gamma, interval):

		self.MODE = MODE

		# parameters
		self.alpha = alpha
		self.gamma = gamma

		# rewards
		self.r_on   = -1 * interval/100
		self.r_off  = 0 * interval
		self.r_comp = 10000

		# epsilon-greedy
		self.epsilon = 0.2
		self.normal  = 0.1

		# log
		self.logf = open('q_logs.txt', 'a')


	def initializeQ(self):
		# only snow accumulation Q[action][Slevel]
		if(self.MODE==0):
			Qtable = np.zeros((2, 10))

		# snow accumulation and temperature Q[action][Slevel, Tlevel]
		elif(self.MODE==1):
			Qtable = np.zeros(((2, 10, 12)))

		# invalid MODE
		else:
			print('invalid MODE error!')
			sys.exit()

		comp = False
		return Qtable, comp


	def Tlevel(self, temp):
		if( temp < -8 ):	Tlevel = 0
		elif( temp < -6 ):	Tlevel = 1
		elif( temp < -4 ):	Tlevel = 2
		elif( temp < -3 ):	Tlevel = 3
		elif( temp < -2 ):	Tlevel = 4
		elif( temp < -1 ):	Tlevel = 5
		elif( temp < 0 ):	Tlevel = 6
		elif( temp < 1 ):	Tlevel = 7
		elif( temp < 2 ):	Tlevel = 8
		elif( temp < 3 ):	Tlevel = 9
		elif( temp < 4 ):	Tlevel = 10
		else:				Tlevel = 11
		return Tlevel


	def Slevel(self, snow):
		if( snow < 0 ):
			print('invalid value of snow accumulation')
			sys.exit()
		elif( snow < 0.01 ):	Slevel = 0
		elif( snow < 0.02 ):	Slevel = 1
		elif( snow < 0.03 ):	Slevel = 2
		elif( snow < 0.04 ):	Slevel = 3
		elif( snow < 0.05 ):	Slevel = 4
		elif( snow < 0.07 ):	Slevel = 5
		elif( snow < 0.1 ) :	Slevel = 6
		elif( snow < 0.15 ):	Slevel = 7
		elif( snow < 0.2 ) :	Slevel = 8
		else:					Slevel = 9
		return Slevel


	def nextMax0(self, Q, heater, eachSLv):
		if(heater==1):
			if( Q[1][eachSLv] > Q[0][eachSLv] ):
				return Q[1][eachSLv]
			else:
				return Q[0][eachSLv]
		else:
			if( Q[1][eachSLv] > Q[0][eachSLv] ):
				return Q[1][eachSLv]
			else:
				return Q[0][eachSLv]


	def nextMax1(self, Q, heater, eachSLv, nextTLv):
		if(heater==1):
			if( Q[1][eachSLv][nextTLv] > Q[0][eachSLv][nextTLv] ):
				return Q[1][eachSLv][nextTLv]
			else:
				return Q[0][eachSLv][nextTLv]
		else:
			if( Q[1][eachSLv][nextTLv] > Q[0][eachSLv][nextTLv] ):
				return Q[1][eachSLv][nextTLv]
			else:
				return Q[0][eachSLv][nextTLv]


	def update_Q0(self, Q, comp, heater, Slevel, onSLv, offSLv):
		reward = 0
		if(Slevel<=4 and comp==False):
			reward = self.r_comp
			comp = True
		else:
			comp = False

		if(heater==1):
			reward += self.r_on
			nextMax = self.nextMax0(Q, 1, onSLv)
			Q[1][Slevel] = (1-self.alpha) * Q[1][Slevel] \
							+ self.alpha * (reward + self.gamma*nextMax)
		else:
			nextMax = self.nextMax0(Q, 0, offSLv)
			Q[0][Slevel] = (1-self.alpha) * Q[0][Slevel] \
							+ self.alpha * (reward + self.gamma*nextMax)
		return Q, comp


	def update_Q1(self, Q, comp, heater, Slevel, Tlevel, onSLv, offSLv, nextTLv):
		reward = 0
		if(Slevel<=4 and comp==False):
			reward = self.r_comp
			comp = True
		else:
			comp = False

		if(heater==1):
			reward += self.r_on
			nextMax = self.nextMax1(Q, 1, onSLv, nextTLv)
			Q[1][Slevel][Tlevel] = (1-self.alpha) * Q[1][Slevel][Tlevel] \
								+ self.alpha * (reward + self.gamma*nextMax)
		else:
			nextMax = self.nextMax1(Q, 0, offSLv, nextTLv)
			Q[0][Slevel][Tlevel] = (1-self.alpha) * Q[0][Slevel][Tlevel] \
								+ self.alpha * (reward + self.gamma*nextMax)
		return Q, comp



	def select_act0(self, Q, Slevel, onSLv, offSLv):
		rand = random.random()
		if( rand > (self.epsilon+self.normal) ):
			if( Q[1][onSLv] > Q[0][offSLv] ):
				act = 1
			elif( Q[1][onSLv] == Q[0][offSLv] ):
				print('same next Q -> random choice')
				self.logf.write('rnadom, ')
				act = random.choice(action)
			else:
				act = 0
		elif( rand <= (self.epsilon+self.normal) and rand > self.epsilon ):
			if(Slevel<=4):
				act = 0
			else:
				act = 1	
		else:
			print('random choice (epsilon-greedy)')
			self.logf.write('rnadom, ')
			act = random.choice(action)
		self.logf.write(str(act)+'\n')
		return act


	def select_act1(self, Q, Slevel, onSLv, offSLv, nextTLv):
		rand = random.random()
		if( rand > (self.epsilon+self.normal) ):
			if( Q[1][onSLv][nextTLv] > Q[0][offSLv][nextTLv] ):
				act = 1
			elif( Q[1][onSLv][nextTLv] == Q[0][offSLv][nextTLv] ):
				print('same next Q -> random choice')
				self.logf.write('rnadom, ')
				act = random.choice(action)
			else:
				act = 0
		elif( rand <= (self.epsilon+self.normal) and rand > self.epsilon ):
			if(Slevel<=4):
				act = 0
			else:
				act = 1	
		else:
			print('random choice')
			self.logf.write('rnadom, ')
			act = random.choice(action)
		return act
