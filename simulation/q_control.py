#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
from datetime import datetime
import argparse
import sys
import numpy as np
import random
import time


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
		self.r_on   = -1 * interval/10
		self.r_comp = 1000

		# epsilon-greedy
		self.epsilon = 0.2
		self.normal  = 0.0

		# log
		self.logf = open('q_logs_'+MODE+'.txt', 'a')


	def initializeQ(self):
		# accumulation Q[action][Slevel]
		if( self.MODE == 'S' ):
			Qtable = np.zeros((2, 7))

		# snow accumulation and temperature Q[action][Slevel][Tlevel]
		elif( self.MODE == 'ST' ):
			Qtable = np.zeros(((2, 7, 6)))

		# precipitation Q[action][Plevel]
		elif( self.MODE == 'P' ):
			Qtable = np.zeros((2, 5))

		# temperature and Precipitation Q[action][Tlevel][Plevel]
		elif( self.MODE == 'TP' ):
			Qtable = np.zeros(((2, 6, 5)))

		# snow accumulation and precipitation Q[action][Slevel][Plevel]
		elif( self.MODE=='SP' ):
			Qtable = np.zeros(((2, 7, 5)))

		# invalid MODE
		else:
			print('invalid MODE error')
			sys.exit()

		comp = False
		return Qtable, comp


	def Tlevel(self, temp):
		if( temp < -6 ):	Tlevel = 0
		elif( temp < -3 ):	Tlevel = 1
		elif( temp < 0 ):	Tlevel = 2
		elif( temp < 3 ):	Tlevel = 3
		elif( temp < 6 ):	Tlevel = 4
		else:				Tlevel = 5
		return Tlevel


	def Slevel(self, snow):
		if( snow < 0 ):
			print('invalid value of snow accumulation')
			sys.exit()
		elif( snow < 0.03 ):	Slevel = 0
		elif( snow < 0.05 ):	Slevel = 1
		elif( snow < 0.1 ) :	Slevel = 2
		elif( snow < 0.2 ) :	Slevel = 3
		elif( snow < 0.3 ) :	Slevel = 4
		elif( snow < 0.5 ) :	Slevel = 5
		else:					Slevel = 6
		return Slevel

	def Plevel(self, pre):
		if( pre < 0 ):
			print('invalid value of snow accumulation')
			sys.exit()
		elif( pre == 0 ):	Plevel = 0
		elif( pre < 0.01 ):	Plevel = 1
		elif( pre < 0.02 ):	Plevel = 2
		elif( pre < 0.03 ):	Plevel = 3
		else:				Plevel = 4
		return Plevel


	def nextMax_S(self, Q, eachSLv):
		if( Q[1][eachSLv] > Q[0][eachSLv] ):
			return Q[1][eachSLv]
		else:
			return Q[0][eachSLv]


	def updateQ_S(self, Q, comp, heater, Slevel, onSLv, offSLv):
		reward = 0
		if(Slevel<=2):
			if(comp==False):
				reward += self.r_comp
			comp = True
		else:
			comp = False

		if(heater==1):
			reward += self.r_on
			nextMax = self.nextMax_S(Q, onSLv)
			Q[1][Slevel] = (1-self.alpha) * Q[1][Slevel] \
							+ self.alpha * (reward + self.gamma*nextMax)
			Q[1][Slevel] = '{:.5f}' .format(Q[1][Slevel])
		else:
			nextMax = self.nextMax_S(Q, offSLv)
			Q[0][Slevel] = (1-self.alpha) * Q[0][Slevel] \
							+ self.alpha * (reward + self.gamma*nextMax)
			Q[0][Slevel] = '{:.5f}' .format(Q[0][Slevel])
		return Q, comp


	def nextMax_ST(self, Q, eachSLv, nextTLv):
		if( Q[1][eachSLv][nextTLv] > Q[0][eachSLv][nextTLv] ):
			return Q[1][eachSLv][nextTLv]
		else:
			return Q[0][eachSLv][nextTLv]


	def updateQ_ST(self, Q, comp, heater, Slevel, Tlevel, onSLv, offSLv, nextTLv):
		reward = 0
		if(Slevel<=2):
			if(comp==False):
				reward += self.r_comp
			comp = True
		else:
			comp = False

		if(heater==1):
			reward += self.r_on
			nextMax = self.nextMax_ST(Q, onSLv, nextTLv)
			Q[1][Slevel][Tlevel] = (1-self.alpha) * Q[1][Slevel][Tlevel] \
								+ self.alpha * (reward + self.gamma*nextMax)
			Q[1][Slevel][Tlevel] = '{:.3f}' .format(Q[1][Slevel][Tlevel])
		else:
			nextMax = self.nextMax_ST(Q, offSLv, nextTLv)
			Q[0][Slevel][Tlevel] = (1-self.alpha) * Q[0][Slevel][Tlevel] \
								+ self.alpha * (reward + self.gamma*nextMax)
			Q[0][Slevel][Tlevel] = '{:.3f}' .format(Q[0][Slevel][Tlevel])
		return Q, comp


	def nextMax_P(self, Q, nextPLv):
		if( Q[1][nextPLv] > Q[0][nextPLv] ):
			return Q[1][nextPLv]
		else:
			return Q[0][nextPLv]


	def updateQ_P(self, Q, comp, heater, Slevel, Plevel, nextPLv):
		reward = 0
		if(Slevel<=2):
			if(comp==False):
				reward += self.r_comp
			comp = True
		else:
			comp = False

		nextMax = self.nextMax_P(Q, nextPLv)
		if(heater==1):
			reward += self.r_on
			Q[1][Plevel] = (1-self.alpha) * Q[1][Plevel]\
							+ self.alpha * (reward + self.gamma*nextMax)
			Q[1][Plevel] = '{:.5f}' .format(Q[1][Plevel])
		else:
			Q[0][Plevel] = (1-self.alpha) * Q[0][Plevel]\
							+ self.alpha * (reward + self.gamma*nextMax)
			Q[0][Plevel] = '{:.5f}' .format(Q[0][Plevel])
		return Q, comp


	def nextMax_TP(self, Q, nextTLv, nextPLv):
		if( Q[1][nextTLv][nextPLv] > Q[0][nextTLv][nextPLv] ):
			return Q[1][nextTLv][nextPLv]
		else:
			return Q[0][nextTLv][nextPLv]


	def updateQ_TP(self, Q, comp, heater, Slevel, Tlevel, Plevel, nextTLv, nextPLv):
		reward = 0
		if(Slevel<=2):
			if(comp==False):
				reward += self.r_comp
			comp = True
		else:
			comp = False

		nextMax = self.nextMax_TP(Q, nextTLv, nextPLv)
		if(heater==1):
			reward += self.r_on
			Q[1][Tlevel][Plevel] = (1-self.alpha) * Q[1][Tlevel][Plevel]\
								+ self.alpha * (reward + self.gamma*nextMax)
			Q[1][Tlevel][Plevel] = '{:.3f}' .format(Q[1][Tlevel][Plevel])
		else:
			Q[0][Tlevel][Plevel] = (1-self.alpha) * Q[0][Tlevel][Plevel]\
								+ self.alpha * (reward + self.gamma*nextMax)
			Q[0][Tlevel][Plevel] = '{:.3f}' .format(Q[0][Tlevel][Plevel])
		return Q, comp


	def nextMax_SP(self, Q, eachSLv, nextPLv):
		if( Q[1][eachSLv][nextPLv] > Q[0][eachSLv][nextPLv] ):
			return Q[1][eachSLv][nextPLv]
		else:
			return Q[0][eachSLv][nextPLv]


	def updateQ_SP(self, Q, comp, heater, Slevel, Plevel, onSLv, offSLv, nextPLv):
		reward = 0
		if(Slevel<=2):
			if(comp==False):
				reward += self.r_comp
			comp = True
		else:
			comp = False

		if(heater==1):
			reward += self.r_on
			nextMax = self.nextMax_SP(Q, onSLv, nextPLv)
			Q[1][Slevel][Plevel] = (1-self.alpha) * Q[1][Slevel][Plevel] \
								+ self.alpha * (reward + self.gamma*nextMax)
			Q[1][Slevel][Plevel] = '{:.3f}' .format(Q[1][Slevel][Plevel])
		else:
			nextMax = self.nextMax_SP(Q, offSLv, nextPLv)
			Q[0][Slevel][Plevel] = (1-self.alpha) * Q[0][Slevel][Plevel] \
								+ self.alpha * (reward + self.gamma*nextMax)
			Q[0][Slevel][Plevel] = '{:.3f}' .format(Q[0][Slevel][Plevel])
		return Q, comp


	def random_act(self, rand, Slevel):
		if( rand <= (self.epsilon+self.normal) and rand > self.epsilon ):
			self.logf.write(',\tnormal, ')
			if(Slevel<=2):		act = 0
			else:				act = 1
		else:
			self.logf.write(',\trandom, ')
			act = random.choice(action)
		return act


	def select_act_S(self, Q, Slevel, onSLv, offSLv):
		rand = random.random()
		if( rand > (self.epsilon+self.normal) ):
			if( Q[1][onSLv] > Q[0][offSLv] ):
				self.logf.write(',\tnormal, ')
				act = 1
			else:
				self.logf.write(',\tnormal, ')
				act = 0
		else:
			act = self.random_act(rand, Slevel)
		self.logf.write('\t'+str(act)+',\t')
		return act


	def select_act_ST(self, Q, Slevel, Tlevel, onSLv, offSLv, nextTLv):
		rand = random.random()
		if( rand > (self.epsilon+self.normal) ):
			if( Q[1][Slevel][Tlevel] > Q[0][Slevel][Tlevel] ):
				self.logf.write(',\tnormal, ')
				act = 1
			else:
				self.logf.write(',\tnormal, ')
				act = 0
		else:
			act = self.random_act(rand, Slevel)
		self.logf.write('\t'+str(act)+',\t')
		return act


	def select_act_P(self, Q, Slevel, Plevel, nextPLv):
		rand = random.random()
		if( rand > (self.epsilon+self.normal) ):
			if( Q[1][Plevel] > Q[0][Plevel] ):
				self.logf.write(',\tnormal, ')
				act = 1
			else:
				self.logf.write(',\tnormal, ')
				act = 0
		else:
			act = self.random_act(rand, Slevel)
		self.logf.write('\t'+str(act)+',\t')
		return act


	def select_act_TP(self, Q, Slevel, Tlevel, Plevel, nextTLv, nextPLv):
		rand = random.random()
		if( rand > (self.epsilon+self.normal) ):
			if( Q[1][Tlevel][Plevel] > Q[0][Tlevel][Plevel] ):
				self.logf.write(',\tnormal, ')
				act = 1
			else:
				self.logf.write(',\tnormal, ')
				act = 0
		else:
			act = self.random_act(rand, Slevel)
		self.logf.write('\t'+str(act)+',\t')
		return act


	def select_act_SP(self, Q, Slevel, Plevel, onSLv, offSLv, nextPLv):
		rand = random.random()
		if( rand > (self.epsilon+self.normal) ):
			if( Q[1][Slevel][Plevel] > Q[0][Slevel][Plevel] ):
				self.logf.write(',\tnormal, ')
				act = 1
			else:
				self.logf.write(',\tnormal, ')
				act = 0
		else:
			act = self.random_act(rand, Slevel)
		self.logf.write('\t'+str(act)+',\t')
		return act



if __name__ == '__main__':
	start = time.time()
	print('start time :', time.ctime())

	action = ['north', 'south', 'east', 'west']
	alpha = 0.1
	gamma = 0.9
	epsilon = 0.2

	# rewards
	r_move = -1
	r_goal = 100

	# initialize Q table
	Q = np.zeros(((4, 4, 4)))


	try:
		for num in range(10000):
			print('----- episode {} -----' .format(num+1))
			state = np.zeros(2, dtype=int)
			nextstate = np.zeros(2, dtype=int)

			while(True):

				# decide action
				if( random.random() > epsilon ):
					if( state[1]<3 ):
						Qmax = Q[state[0]][state[1]][0]
						act = 'north'
					elif( state[1]>0 ):
						Qmax = Q[state[0]][state[1]][1]
						act = 'south'

					if( Q[state[0]][state[1]][1]>Qmax and state[1]>0 ):
						Qmax = Q[state[0]][state[1]][1]
						act = 'south'
					if( Q[state[0]][state[1]][2]>Qmax and state[0]<3 ):
						Qmax = Q[state[0]][state[1]][2]
						act = 'east'
					if( Q[state[0]][state[1]][3]>Qmax and state[0]>0 ):
						Qmax = Q[state[0]][state[1]][3]
						act = 'west'

				else:
					del_list = []
					if( state[1]>=3 ):
						del_list.append(0)
					elif( state[1]<=0 ):
						del_list.append(1)
					if( state[0]>=3 ):
						del_list.append(2)
					elif( state[0]<=0 ):
						del_list.append(3)
					A = np.delete(action, del_list)
					act = random.choice(A)

				if( act=='north' ):
					a = 0
					nextstate[1] += 1
					print(act)
				elif( act=='south' ):
					a = 1
					nextstate[1] -= 1
					print(act)
				elif( act=='east' ):
					a = 2
					nextstate[0] += 1
					print(act)
				else:
					a = 3
					nextstate[0] -= 1
					print(act)
				print(state)

				# update Q table
				reward = r_move
				if( nextstate[0]==3 and nextstate[1]==3 ):
					reward += r_goal
					print('goal!')

				nextMax = max(Q[nextstate[0]][nextstate[1]])
				Q[state[0]][state[1]][a] = (1-alpha) * Q[state[0]][state[1]][a]\
											+ alpha * (reward + gamma*nextMax)

	
				if( act=='north' ):		state[1] += 1
				elif( act=='south' ):	state[1] -= 1
				elif( act=='east' ):	state[0] += 1
				else:					state[0] -= 1

				if( state[0]==3 and state[1]==3 ):
					break

#			time.sleep(1)
#			print(Q)

	finally:
		print(Q)
