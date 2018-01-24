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

	def __init__(self, MODE, Qloop, alpha, gamma, interval):
		self.MODE = MODE

		# parameters
		self.alpha = alpha
		self.gamma = gamma

		# rewards
		self.r_on   = -1 * interval
		self.r_off  = 1 * interval/10
		self.r_comp = 100

		# epsilon-greedy
		self.epsilon = 0.2
		self.normal  = 0.0

		# log
		self.logf = open('Qlogs/q_logs_'+MODE+str(Qloop)+'.txt', 'a')


	def initializeQ(self):
		# snow accumulation Q[action][Srank]
		if( self.MODE == 'S' ):
			Qtable = np.zeros((2, 7))

		# snow accumulation and temperature Q[action][Srank][Trank]
		elif( self.MODE == 'ST' ):
			Qtable = np.zeros(((2, 7, 6)))

		# precipitation Q[action][Prank]
		elif( self.MODE == 'P' ):
			Qtable = np.zeros((2, 5))

		# temperature and Precipitation Q[action][Trank][Prank]
		elif( self.MODE == 'TP' ):
			Qtable = np.zeros(((2, 6, 5)))

		# snow accumulation and precipitation Q[action][Srank][Prank]
		elif( self.MODE=='SP' ):
			Qtable = np.zeros(((2, 7, 5)))

		# continuously operating time Q[action][ONtime]
		elif( self.MODE=='C' ):
			Qtable = np.zeros((2, 5))

		# snow accumulation & continuously operating time
		# Q[action][Srank][Crank]
		elif( self.MODE=='SC' ):
			Qtable = np.zeros(((2, 7, 5)))

		# precipitation and continuously operating time Q[action][Plevel][Clevel]
		elif( self.MODE=='PC' ):
			Qtable = np.zeros(((2, 5, 5)))

		# invalid MODE
		else:
			print('invalid MODE error')
			sys.exit()

		comp = False
		return Qtable, comp


	def Trank(self, temp):
		if( temp < -6 ):	Trank = 0
		elif( temp < -3 ):	Trank = 1
		elif( temp < 0 ):	Trank = 2
		elif( temp < 3 ):	Trank = 3
		elif( temp < 6 ):	Trank = 4
		else:				Trank = 5
		return Trank


	def Srank(self, snow):
		if( snow < 0 ):
			print('invalid value of snow accumulation')
			sys.exit()
		elif( snow < 0.03 ):	Srank = 0
		elif( snow < 0.05 ):	Srank = 1
		elif( snow < 0.1 ) :	Srank = 2
		elif( snow < 0.2 ) :	Srank = 3
		elif( snow < 0.4 ) :	Srank = 4
		elif( snow < 0.7 ) :	Srank = 5
		else:					Srank = 6
		return Srank

	def Prank(self, pre):
		if( pre < 0 ):
			print('invalid value of snow accumulation')
			sys.exit()
		elif( pre == 0 ):	Prank = 0
		elif( pre < 0.01 ):	Prank = 1
		elif( pre < 0.02 ):	Prank = 2
		elif( pre < 0.03 ):	Prank = 3
		else:				Prank = 4
		return Prank

	
	def Crank(self, onT):
		if( onT == 0 ):		Crank = 0
		elif( onT < 10 ):	Crank = 1
		elif( onT < 20 ):	Crank = 2
		elif( onT < 30 ):	Crank = 3
		else:				Crank = 4
		return Crank


	def nextMax_S(self, Q, eachSr):
		if( Q[1][eachSr] > Q[0][eachSr] ):
			return Q[1][eachSr]
		else:
			return Q[0][eachSr]


	def updateQ_S(self, Q, comp, heater, Srank, onSr, offSr):
		reward = 0
		if(Srank<=1):
			if(comp==False):
				reward += self.r_comp
			comp = True
		else:
			comp = False

		if(heater==1):
			reward += self.r_on
			nextMax = self.nextMax_S(Q, onSr)
			Q[1][Srank] = (1-self.alpha) * Q[1][Srank] \
							+ self.alpha * (reward + self.gamma*nextMax)
			Q[1][Srank] = '{:.5f}' .format(Q[1][Srank])
		else:
			reward += self.r_off
			nextMax = self.nextMax_S(Q, offSr)
			Q[0][Srank] = (1-self.alpha) * Q[0][Srank] \
							+ self.alpha * (reward + self.gamma*nextMax)
			Q[0][Srank] = '{:.5f}' .format(Q[0][Srank])
		return Q, comp


	def nextMax_ST(self, Q, eachSr, nextTr):
		if( Q[1][eachSr][nextTr] > Q[0][eachSr][nextTr] ):
			return Q[1][eachSr][nextTr]
		else:
			return Q[0][eachSr][nextTr]


	def updateQ_ST(self, Q, comp, heater, Srank, Trank, onSr, offSr, nextTr):
		reward = 0
		if(Srank<=1):
			if(comp==False):
				reward += self.r_comp
			comp = True
		else:
			comp = False

		if(heater==1):
			reward += self.r_on
			nextMax = self.nextMax_ST(Q, onSr, nextTr)
			Q[1][Srank][Trank] = (1-self.alpha) * Q[1][Srank][Trank] \
								+ self.alpha * (reward + self.gamma*nextMax)
			Q[1][Srank][Trank] = '{:.3f}' .format(Q[1][Srank][Trank])
		else:
			reward += self.r_off
			nextMax = self.nextMax_ST(Q, offSr, nextTr)
			Q[0][Srank][Trank] = (1-self.alpha) * Q[0][Srank][Trank] \
								+ self.alpha * (reward + self.gamma*nextMax)
			Q[0][Srank][Trank] = '{:.3f}' .format(Q[0][Srank][Trank])
		return Q, comp


	def nextMax_P(self, Q, nextPr):
		if( Q[1][nextPr] > Q[0][nextPr] ):
			return Q[1][nextPr]
		else:
			return Q[0][nextPr]


	def updateQ_P(self, Q, comp, heater, Srank, Prank, nextPr):
		reward = 0
		if(Srank<=1):
			if(comp==False):
				reward += self.r_comp
			comp = True
		else:
			comp = False

		nextMax = self.nextMax_P(Q, nextPr)
		if(heater==1):
			reward += self.r_on
			Q[1][Prank] = (1-self.alpha) * Q[1][Prank]\
							+ self.alpha * (reward + self.gamma*nextMax)
			Q[1][Prank] = '{:.5f}' .format(Q[1][Prank])
		else:
			reward += self.r_off
			Q[0][Prank] = (1-self.alpha) * Q[0][Prank]\
							+ self.alpha * (reward + self.gamma*nextMax)
			Q[0][Prank] = '{:.5f}' .format(Q[0][Prank])
		return Q, comp


	def nextMax_C(self, Q, nextCr):
		if( Q[1][nextCr] > Q[0][nextCr] ):
			return Q[1][nextCr]
		else:
			return Q[0][nextCr]


	def updateQ_C(self, Q, comp, heater, Srank, Crank, onCr, offCr):
		reward = 0
		if(Srank<=1):
			if(comp==False):
				reward += self.r_comp
			comp = True
		else:
			comp = False

		if(heater==1):
			reward += self.r_on
			nextMax = self.nextMax_C(Q, onCr)
			Q[1][Crank] = (1-self.alpha) * Q[1][Crank] \
							+ self.alpha * (reward + self.gamma*nextMax)
			Q[1][Crank] = '{:.5f}' .format(Q[1][Crank])
		else:
			reward += self.r_off
			nextMax = self.nextMax_C(Q, offCr)
			Q[0][Crank] = (1-self.alpha) * Q[0][Crank] \
							+ self.alpha * (reward + self.gamma*nextMax)
			Q[0][Crank] = '{:.5f}' .format(Q[0][Crank])
		return Q, comp


	def nextMax_TP(self, Q, nextTr, nextPr):
		if( Q[1][nextTr][nextPr] > Q[0][nextTr][nextPr] ):
			return Q[1][nextTr][nextPr]
		else:
			return Q[0][nextTr][nextPr]


	def updateQ_TP(self, Q, comp, heater, Srank, Trank, Prank, nextTr, nextPr):
		reward = 0
		if(Srank<=1):
			if(comp==False):
				reward += self.r_comp
			comp = True
		else:
			comp = False

		nextMax = self.nextMax_TP(Q, nextTr, nextPr)
		if(heater==1):
			reward += self.r_on
			Q[1][Trank][Prank] = (1-self.alpha) * Q[1][Trank][Prank]\
								+ self.alpha * (reward + self.gamma*nextMax)
			Q[1][Trank][Prank] = '{:.3f}' .format(Q[1][Trank][Prank])
		else:
			reward += self.r_off
			Q[0][Trank][Prank] = (1-self.alpha) * Q[0][Trank][Prank]\
								+ self.alpha * (reward + self.gamma*nextMax)
			Q[0][Trank][Prank] = '{:.3f}' .format(Q[0][Trank][Prank])
		return Q, comp


	def nextMax_SP(self, Q, eachSr, nextPr):
		if( Q[1][eachSr][nextPr] > Q[0][eachSr][nextPr] ):
			return Q[1][eachSr][nextPr]
		else:
			return Q[0][eachSr][nextPr]


	def updateQ_SP(self, Q, comp, heater, Srank, Prank, onSr, offSr, nextPr):
		reward = 0
		if(Srank<=1):
			if(comp==False):
				reward += self.r_comp
			comp = True
		else:
			comp = False

		if(heater==1):
			reward += self.r_on
			nextMax = self.nextMax_SP(Q, onSr, nextPr)
			Q[1][Srank][Prank] = (1-self.alpha) * Q[1][Srank][Prank] \
								+ self.alpha * (reward + self.gamma*nextMax)
			Q[1][Srank][Prank] = '{:.3f}' .format(Q[1][Srank][Prank])
		else:
			reward += self.r_off
			nextMax = self.nextMax_SP(Q, offSr, nextPr)
			Q[0][Srank][Prank] = (1-self.alpha) * Q[0][Srank][Prank] \
								+ self.alpha * (reward + self.gamma*nextMax)
			Q[0][Srank][Prank] = '{:.3f}' .format(Q[0][Srank][Prank])
		return Q, comp


	def nextMax_SC(self, Q, eachSr, eachCr):
		if( Q[1][eachSr][eachCr] > Q[0][eachSr][eachCr] ):
			return Q[1][eachSr][eachCr]
		else:
			return Q[0][eachSr][eachCr]


	def updateQ_SC(self, Q, comp, heater, Srank, Crank,\
					onSr, offSr, onCr, offCr):
		reward = 0
		if(Srank<=1):
			if(comp==False):
				reward += self.r_comp
			else:
				comp = True
		else:
			comp = False

		if(heater==1):
			reward += self.r_on
			nextMax = self.nextMax_SC(Q, onSr, onCr)
			Q[1][Srank][Crank] = (1-self.alpha) * Q[1][Srank][Crank]\
									+ self.alpha * (reward + self.gamma*nextMax)
			Q[1][Srank][Crank] = '{:.3f}' .format(Q[1][Srank][Crank])
		else:
			reward += self.r_off
			nextMax = self.nextMax_SC(Q, offSr, offCr)
			Q[0][Srank][Crank] = (1-self.alpha) * Q[0][Srank][Crank]\
									+ self.alpha * (reward + self.gamma*nextMax)
			Q[0][Srank][Crank] = '{:.3f}' .format(Q[0][Srank][Crank])
		return Q, comp


	def nextMax_PC(self, Q, nextPr, nextCr):
		if( Q[1][nextPr][nextCr] > Q[0][nextPr][nextCr] ):
			return Q[1][nextPr][nextCr]
		else:
			return Q[0][nextPr][nextCr]


	def updateQ_PC(self, Q, comp, heater, Srank, Prank, Crank, nextPr, onCr, offCr):
		reward = 0
		if(Srank<=1):
			if(comp==False):
				reward += self.r_comp
			comp = True
		else:
			comp = False

		if(heater==1):
			nextMax = self.nextMax_PC(Q, nextPr, onCr)
			reward += self.r_on
			Q[1][Prank][Crank] = (1-self.alpha) * Q[1][Prank][Crank]\
								+ self.alpha * (reward + self.gamma*nextMax)
			Q[1][Prank][Crank] = '{:.3f}' .format(Q[1][Prank][Crank])
		else:
			reward += self.r_off
			nextMax = self.nextMax_PC(Q, nextPr, offCr)
			Q[0][Prank][Crank] = (1-self.alpha) * Q[0][Prank][Crank]\
								+ self.alpha * (reward + self.gamma*nextMax)
			Q[0][Prank][Crank] = '{:.3f}' .format(Q[0][Prank][Crank])
		return Q, comp


	def random_act(self, rand, Srank):
		if( rand <= (self.epsilon+self.normal) and rand > self.epsilon ):
			self.logf.write(',\tnormal,')
			if(Srank<=1):		act = 0
			else:				act = 1
		else:
			self.logf.write(',\trandom,')
			act = random.choice(action)
		return act


	def select_act_S(self, Q, Srank):
		rand = random.random()
		if( rand > (self.epsilon+self.normal) ):
			if( Q[1][Srank] > Q[0][Srank] ):
				self.logf.write(',\tnormal,')
				act = 1
			else:
				self.logf.write(',\tnormal,')
				act = 0
		else:
			act = self.random_act(rand, Srank)
		self.logf.write('\t'+str(act)+',\t')
		return act


	def select_act_ST(self, Q, Srank, Trank):
		rand = random.random()
		if( rand > (self.epsilon+self.normal) ):
			if( Q[1][Srank][Trank] > Q[0][Srank][Trank] ):
				self.logf.write(',\tnormal,')
				act = 1
			else:
				self.logf.write(',\tnormal,')
				act = 0
		else:
			act = self.random_act(rand, Srank)
		self.logf.write('\t'+str(act)+',\t')
		return act


	def select_act_P(self, Q, Srank, Prank):
		rand = random.random()
		if( rand > (self.epsilon+self.normal) ):
			if( Q[1][Prank] > Q[0][Prank] ):
				self.logf.write(',\tnormal,')
				act = 1
			else:
				self.logf.write(',\tnormal,')
				act = 0
		else:
			act = self.random_act(rand, Srank)
		self.logf.write('\t'+str(act)+',\t')
		return act


	def select_act_TP(self, Q, Srank, Trank, Prank):
		rand = random.random()
		if( rand > (self.epsilon+self.normal) ):
			if( Q[1][Trank][Prank] > Q[0][Trank][Prank] ):
				self.logf.write(',\tnormal,')
				act = 1
			else:
				self.logf.write(',\tnormal,')
				act = 0
		else:
			act = self.random_act(rand, Srank)
		self.logf.write('\t'+str(act)+',\t')
		return act


	def select_act_SP(self, Q, Srank, Prank):
		rand = random.random()
		if( rand > (self.epsilon+self.normal) ):
			if( Q[1][Srank][Prank] > Q[0][Srank][Prank] ):
				self.logf.write(',\tnormal,')
				act = 1
			else:
				self.logf.write(',\tnormal,')
				act = 0
		else:
			act = self.random_act(rand, Srank)
		self.logf.write('\t'+str(act)+',\t')
		return act


	def select_act_C(self, Q, Srank, Crank, onCr, offCr):
		rand = random.random()
		if( rand > (self.epsilon+self.normal) ):
			if( Q[1][Crank] > Q[0][Crank] ):
				self.logf.write(',\tnormal,')
				act = 1
			else:
				self.logf.write(',\tnormal,')
				act = 0
		else:
			act = self.random_act(rand, Srank)
		self.logf.write('\t'+str(act)+',\t')
		return act


	def select_act_SC(self, Q, Srank, Crank):
		rand = random.random()
		if( rand > (self.epsilon+self.normal) ):
			if( Q[1][Srank][Crank] > Q[0][Srank][Crank] ):
				self.logf.write(',\tnormal,')
				act = 1
			else:
				self.logf.write(',\tnormal,')
				act = 0
		else:
			act = self.random_act(rand, Srank)
		self.logf.write('\t'+str(act)+',\t')
		return act


	def select_act_PC(self, Q, Srank, Prank, Crank):
		rand = random.random()
		if( rand > (self.epsilon+self.normal) ):
			if( Q[1][Prank][Crank] > Q[0][Prank][Crank] ):
				self.logf.write(',\tnormal,')
				act = 1
			else:
				self.logf.write(',\tnormal,')
				act = 0
		else:
			act = self.random_act(rand, Srank)
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
