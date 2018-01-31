#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
import datetime
import argparse
import sys
import numpy as np
import random
import time

# time list of start long snowing
ssList = open('ssList.txt', 'r').readlines()
start_list = []
for ss in ssList:
	st = ss.split('\n')[0]
	start_list.append(datetime.datetime.strptime(st, '%Y-%m-%d %H:%M'))



class control:

	# control from snowfall sensor
	def judge_p(self, pre):
		if(pre > 0):		heater = 1
		else:				heater = 0
		return heater


	# control from snow accumulation
	def judge_s(self, snow):
		if(snow >= 0.2):	heater = 1
		else:				heater = 0
		return heater


	# always off
	def off(self):
		heater = 0
		return heater


	# always on
	def on(self):
		heater = 1
		return heater


	# melt after the snowing stops
	def original1(self, pre, snow):
		if( snow < 0.2 ):
			heater = 0
		elif( pre > 0.0 ):
			heater = 0
		else:
			heater = 1
		return heater


	# keep the road temperature above 10
	def original2_0(self, TS):
		if( TS < 10 ):
			heater = 1
		else:
			heater = 0
		return heater


	# keep the road temperature above 0  version.2
	def original2_1(self, TS, snow):
		if( (TS < 0) or (snow >= 0.2) ):
			heater = 1
		else:
			heater = 0
		return heater


	# keep the road temperature above 10  version.2
	def original2_2(self, TS, snow):
		if( (TS < 10) or (snow >= 0.2) ):
			heater = 1
		else:
			heater = 0
		return heater


	# start warming up the road before snow falls continuously
	def original3_0(self, date, snow, TS):
		for start in start_list:
			back = start - datetime.timedelta(minutes=20)
			if( date>=back and date<=start and TS<60 ):
				heater = 1
				break
			else:
				heater = 0
		if( snow >= 0.2 ):
			heater = 1
		return heater


	# start warming up the road before snow falls continuously version.2
	def original3_1(self, date, snow, TS, pre):
		for start in start_list:
			back = start - datetime.timedelta(minutes=20)
			if( date>=back and date<=start and TS<60 ):
				heater = 1
				break
			else:
				heater = 0
		if( pre > 0.0 and TS<60 ):
			heater = 1
		if( snow >= 0.2 ):
			heater = 1
		return heater


	# switch on/off at fixed time while snowing and melt all after snow stops
	def original4_0(self, onT, offT, pre, heater, snow):
		if( pre > 0.0 ):
			if( offT >= 10.0 ):
				heater = 1
			elif( onT >= 10.0 ):
				heater = 0
			elif( onT==0 and offT==0 ):
				heater = 1
		else:
			if( snow >= 0.2 ):
				heater = 1
			else:
				heater = 0
		return heater


	# switch on/off at fixed time while snowing and melt all after snow stops
	# version.2
	def original4_1(self, onT, offT, pre, heater, snow):
		if( pre > 0.0 ):
			if( offT >= 10.0 ):
				heater = 1
			elif( onT >= 20.0 ):
				heater = 0
			elif( onT==0 and offT==0 ):
				heater = 1
		else:
			if( snow >= 0.2 ):
				heater = 1
			else:
				heater = 0
		return heater


	# turn off a little earlier
	def original5_0(self, snow, previous_snow):
		if( snow >= 0.2 ):
			heater = 1
			if( (snow<0.3) and (previous_snow>snow) ):
				heater = 0
		else:
			heater = 0
		return heater


	# turn off a little more earlier
	def original5_1(self, snow, previous_snow):
		if( snow >= 0.2 ):
			heater = 1
			if( (snow<0.5) and (previous_snow>snow) ):
				heater = 0
		else:
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
		self.r_on   = -1 * interval /100.0
		self.r_comp = 1000
		self.r_much = -1 * interval /10.0

		# epsilon-greedy
		self.eps = 0.2
		self.normal = 0.0

		# log
		self.slogf = open('Qlogs/state_log_'+MODE+'.csv', 'a')
		self.qlogf = open('Qlogs/q_log_'+MODE+'.csv', 'a')


	def initializeQ(self):
		# snow accumulation Q[action][Srank]
		if( self.MODE == 'S' ):
			Qtable = np.empty((2, 7))

		# snow accumulation and temperature Q[action][Srank][Trank]
		elif( self.MODE == 'ST' ):
			Qtable = np.empty((2, 7, 6))

		# precipitation Q[action][Prank]
		elif( self.MODE == 'P' ):
			Qtable = np.empty((2, 5))

		# temperature and Precipitation Q[action][Trank][Prank]
		elif( self.MODE == 'TP' ):
			Qtable = np.empty((2, 6, 5))

		# snow accumulation and precipitation Q[action][Srank][Prank]
		elif( self.MODE=='SP' ):
			Qtable = np.empty((2, 7, 5))

		# continuously operating time Q[action][Crank]
		elif( self.MODE=='C' ):
			Qtable = np.empty((2, 5))

		# snow accumulation & continuously operating time
		# Q[action][Srank][Crank]
		elif( self.MODE=='SC' ):
			Qtable = np.empty((2, 7, 5))

		# precipitation and continuously operating time Q[action][Prank][Crank]
		elif( self.MODE=='PC' ):
			Qtable = np.empty((2, 5, 5))

		# snow accumulate, precipitation and continuously operating time
		# Q[action][Srank][Prank][Crank]
		elif( self.MODE=='SPC' ):
			Qtable = np.empty((2, 7, 5, 5))

		# invalid MODE
		else:
			print('invalid MODE error')
			sys.exit()

		Qtable = np.round(Qtable, 3)
		self.qlogf.write(str(Qtable))
		self.qlogf.close()

		comp = True
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
		standard = 0.2
		if( snow < 0 ):
			print('invalid value of snow accumulation')
			sys.exit()
		elif( snow < standard ):		Srank = 0
		elif( snow < standard+0.1 ):	Srank = 1
		elif( snow < standard+0.2 ):	Srank = 2
		elif( snow < standard+0.3 ):	Srank = 3
		elif( snow < standard+0.5 ):	Srank = 4
		elif( snow < standard+1.0 ):	Srank = 5
		else:							Srank = 6
		return Srank


	def Prank(self, pre):
		if( pre < 0 ):
			print('invalid value of precipitation')
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
		if(Srank==0):
			if(comp==False):
				reward += self.r_comp
			comp = True
		else:
			if(Srank>=4 and heater==0):
				reward += self.r_much
			comp = False

		if(heater==1):
			reward += self.r_on
			nextMax = self.nextMax_S(Q, onSr)
			Q[1][Srank] = (1-self.alpha) * Q[1][Srank] \
							+ self.alpha * (reward + self.gamma*nextMax)
			Q[1][Srank] = "%.3f" % Q[1][Srank]
		else:
			nextMax = self.nextMax_S(Q, offSr)
			Q[0][Srank] = (1-self.alpha) * Q[0][Srank] \
							+ self.alpha * (reward + self.gamma*nextMax)
			Q[0][Srank] = "%.3f" % Q[0][Srank]

		self.slogf.write('{:>5}, ' .format(reward))
		self.slogf.close()
		return Q, comp


	def nextMax_ST(self, Q, eachSr, nextTr):
		if( Q[1][eachSr][nextTr] > Q[0][eachSr][nextTr] ):
			return Q[1][eachSr][nextTr]
		else:
			return Q[0][eachSr][nextTr]


	def updateQ_ST(self, Q, comp, heater, Srank, Trank, onSr, offSr, nextTr):
		reward = 0
		if(Srank==0):
			if(comp==False):
				reward += self.r_comp
			comp = True
		else:
			if(Srank>=4 and heater==0):
				reward += self.r_much
			comp = False

		if(heater==1):
			reward += self.r_on
			nextMax = self.nextMax_ST(Q, onSr, nextTr)
			Q[1][Srank][Trank] = (1-self.alpha) * Q[1][Srank][Trank] \
								+ self.alpha * (reward + self.gamma*nextMax)
			Q[1][Srank][Trank] = '{:.3f}' .format(Q[1][Srank][Trank])
		else:
			nextMax = self.nextMax_ST(Q, offSr, nextTr)
			Q[0][Srank][Trank] = (1-self.alpha) * Q[0][Srank][Trank] \
								+ self.alpha * (reward + self.gamma*nextMax)
			Q[0][Srank][Trank] = '{:.3f}' .format(Q[0][Srank][Trank])

		self.slogf.write('{:>5}, ' .format(reward))
		self.slogf.close()
		return Q, comp


	def nextMax_P(self, Q, nextPr):
		if( Q[1][nextPr] > Q[0][nextPr] ):
			return Q[1][nextPr]
		else:
			return Q[0][nextPr]


	def updateQ_P(self, Q, comp, heater, Srank, Prank, nextPr):
		reward = 0
		if(Srank==0):
			if(comp==False):
				reward += self.r_comp
			comp = True
		else:
			if(Srank>=4 and heater==0):
				reward += self.r_much
			comp = False

		nextMax = self.nextMax_P(Q, nextPr)
		if(heater==1):
			reward += self.r_on
			Q[1][Prank] = (1-self.alpha) * Q[1][Prank]\
							+ self.alpha * (reward + self.gamma*nextMax)
			Q[1][Prank] = '{:.5f}' .format(Q[1][Prank])
		else:
			Q[0][Prank] = (1-self.alpha) * Q[0][Prank]\
							+ self.alpha * (reward + self.gamma*nextMax)
			Q[0][Prank] = '{:.5f}' .format(Q[0][Prank])

		self.slogf.write('{:>5}, ' .format(reward))
		self.slogf.close()
		return Q, comp


	def nextMax_C(self, Q, nextCr):
		if( Q[1][nextCr] > Q[0][nextCr] ):
			return Q[1][nextCr]
		else:
			return Q[0][nextCr]


	def updateQ_C(self, Q, comp, heater, Srank, Crank, onCr, offCr):
		reward = 0
		if(Srank==0):
			if(comp==False):
				reward += self.r_comp
			comp = True
		else:
			if(Srank>=4 and heater==0):
				reward += self.r_much
			comp = False

		if(heater==1):
			reward += self.r_on
			nextMax = self.nextMax_C(Q, onCr)
			Q[1][Crank] = (1-self.alpha) * Q[1][Crank] \
							+ self.alpha * (reward + self.gamma*nextMax)
			Q[1][Crank] = '{:.5f}' .format(Q[1][Crank])
		else:
			nextMax = self.nextMax_C(Q, offCr)
			Q[0][Crank] = (1-self.alpha) * Q[0][Crank] \
							+ self.alpha * (reward + self.gamma*nextMax)
			Q[0][Crank] = '{:.5f}' .format(Q[0][Crank])

		self.slogf.write('{:>5}, ' .format(reward))
		self.slogf.close()
		return Q, comp


	def nextMax_TP(self, Q, nextTr, nextPr):
		if( Q[1][nextTr][nextPr] > Q[0][nextTr][nextPr] ):
			return Q[1][nextTr][nextPr]
		else:
			return Q[0][nextTr][nextPr]


	def updateQ_TP(self, Q, comp, heater, Srank, Trank, Prank, nextTr, nextPr):
		reward = 0
		if(Srank==0):
			if(comp==False):
				reward += self.r_comp
			comp = True
		else:
			if(Srank>=4 and heater==0):
				reward += self.r_much
			comp = False

		nextMax = self.nextMax_TP(Q, nextTr, nextPr)
		if(heater==1):
			reward += self.r_on
			Q[1][Trank][Prank] = (1-self.alpha) * Q[1][Trank][Prank]\
								+ self.alpha * (reward + self.gamma*nextMax)
			Q[1][Trank][Prank] = '{:.3f}' .format(Q[1][Trank][Prank])
		else:
			Q[0][Trank][Prank] = (1-self.alpha) * Q[0][Trank][Prank]\
								+ self.alpha * (reward + self.gamma*nextMax)
			Q[0][Trank][Prank] = '{:.3f}' .format(Q[0][Trank][Prank])

		self.slogf.write('{:>5}, ' .format(reward))
		self.slogf.close()
		return Q, comp


	def nextMax_SP(self, Q, eachSr, nextPr):
		if( Q[1][eachSr][nextPr] > Q[0][eachSr][nextPr] ):
			return Q[1][eachSr][nextPr]
		else:
			return Q[0][eachSr][nextPr]


	def updateQ_SP(self, Q, comp, heater, Srank, Prank, onSr, offSr, nextPr):
		reward = 0
		if(Srank==0):
			if(comp==False):
				reward += self.r_comp
			comp = True
		else:
			if(Srank>=4 and heater==0):
				reward += self.r_much
			comp = False

		if(heater==1):
			reward += self.r_on
			nextMax = self.nextMax_SP(Q, onSr, nextPr)
			Q[1][Srank][Prank] = (1-self.alpha) * Q[1][Srank][Prank] \
								+ self.alpha * (reward + self.gamma*nextMax)
			Q[1][Srank][Prank] = '{:.3f}' .format(Q[1][Srank][Prank])
		else:
			nextMax = self.nextMax_SP(Q, offSr, nextPr)
			Q[0][Srank][Prank] = (1-self.alpha) * Q[0][Srank][Prank] \
								+ self.alpha * (reward + self.gamma*nextMax)
			Q[0][Srank][Prank] = '{:.3f}' .format(Q[0][Srank][Prank])

		self.slogf.write('{:>5}, ' .format(reward))
		self.slogf.close()
		return Q, comp


	def nextMax_SC(self, Q, eachSr, eachCr):
		if( Q[1][eachSr][eachCr] > Q[0][eachSr][eachCr] ):
			return Q[1][eachSr][eachCr]
		else:
			return Q[0][eachSr][eachCr]


	def updateQ_SC(self, Q, comp, heater, Srank, Crank,\
					onSr, offSr, onCr, offCr):
		reward = 0
		if(Srank==0):
			if(comp==False):
				reward += self.r_comp
			else:
				comp = True
		else:
			if(Srank>=4 and heater==0):
				reward += self.r_much
			comp = False

		if(heater==1):
			reward += self.r_on
			nextMax = self.nextMax_SC(Q, onSr, onCr)
			Q[1][Srank][Crank] = (1-self.alpha) * Q[1][Srank][Crank]\
									+ self.alpha * (reward + self.gamma*nextMax)
			Q[1][Srank][Crank] = '{:.3f}' .format(Q[1][Srank][Crank])
		else:
			nextMax = self.nextMax_SC(Q, offSr, offCr)
			Q[0][Srank][Crank] = (1-self.alpha) * Q[0][Srank][Crank]\
									+ self.alpha * (reward + self.gamma*nextMax)
			Q[0][Srank][Crank] = '{:.3f}' .format(Q[0][Srank][Crank])

		self.slogf.write('{:>5}, ' .format(reward))
		self.slogf.close()
		return Q, comp


	def nextMax_PC(self, Q, nextPr, nextCr):
		if( Q[1][nextPr][nextCr] > Q[0][nextPr][nextCr] ):
			return Q[1][nextPr][nextCr]
		else:
			return Q[0][nextPr][nextCr]


	def updateQ_PC(self, Q, comp, heater, Srank, Prank, Crank, nextPr, onCr, offCr):
		reward = 0
		if(Srank==0):
			if(comp==False):
				reward += self.r_comp
			comp = True
		else:
			if(Srank>=4 and heater==0):
				reward += self.r_much
			comp = False

		if(heater==1):
			nextMax = self.nextMax_PC(Q, nextPr, onCr)
			reward += self.r_on
			Q[1][Prank][Crank] = (1-self.alpha) * Q[1][Prank][Crank]\
								+ self.alpha * (reward + self.gamma*nextMax)
			Q[1][Prank][Crank] = '{:.3f}' .format(Q[1][Prank][Crank])
		else:
			nextMax = self.nextMax_PC(Q, nextPr, offCr)
			Q[0][Prank][Crank] = (1-self.alpha) * Q[0][Prank][Crank]\
								+ self.alpha * (reward + self.gamma*nextMax)
			Q[0][Prank][Crank] = '{:.3f}' .format(Q[0][Prank][Crank])

		self.slogf.write('{:>5}, ' .format(reward))
		self.slogf.close()
		return Q, comp


	def nextMax_SPC(self, Q, nextSr, nextPr, nextCr):
		if( Q[1][nextSr][nextPr][nextCr] > Q[0][nextSr][nextPr][nextCr] ):
			return Q[1][nextSr][nextPr][nextCr]
		else:
			return Q[0][nextSr][nextPr][nextCr]


	def updateQ_SPC(self, Q, comp, heater, Srank, Prank, Crank, onSr, offSr,\
					nextPr, onCr, offCr):
		reward = 0
		if(Srank==0):
			if(comp==False):
				reward += self.r_comp
			comp = True
		else:
			if(Srank>=4 and heater==0):
				reward += self.r_much
			comp = False

		if(heater==1):
			nextMax = self.nextMax_SPC(Q, onSr, nextPr, onCr)
			reward += self.r_on
			Q[1][Srank][Prank][Crank] = (1-self.alpha) * Q[1][Srank][Prank][Crank]\
									+ self.alpha * (reward + self.gamma*nextMax)
			Q[1][Srank][Prank][Crank] = '{:.2f}' .format(Q[1][Srank][Prank][Crank])
		else:
			nextMax = self.nextMax_SPC(Q, offSr, nextPr, offCr)
			Q[0][Srank][Prank][Crank] = (1-self.alpha) * Q[0][Srank][Prank][Crank]\
									+ self.alpha * (reward + self.gamma*nextMax)
			Q[0][Srank][Prank][Crank] = '{:.2f}' .format(Q[0][Srank][Prank][Crank])

		self.slogf.write('{:>5}, ' .format(reward))
		self.slogf.close()
		return Q, comp


	def random_act(self, rand, Srank):
		if( rand <= (self.eps+self.normal) and rand > self.eps ):
			self.slogf.write(', normal, ')
			if(Srank==0):		act = 0
			else:				act = 1
		else:
			self.slogf.write(', random, ')
			act = random.choice(action)
		return act


	def select_act_S(self, Q, Srank):
		rand = random.random()
		if( rand > (self.eps+self.normal) ):
			if( Q[1][Srank] > Q[0][Srank] ):
				self.slogf.write(', normal, ')
				act = 1
			else:
				self.slogf.write(', normal, ')
				act = 0
		else:
			act = self.random_act(rand, Srank)
		self.slogf.write(str(act)+', ')
		self.slogf.close()
		return act


	def select_act_ST(self, Q, Srank, Trank):
		rand = random.random()
		if( rand > (self.eps+self.normal) ):
			if( Q[1][Srank][Trank] > Q[0][Srank][Trank] ):
				self.slogf.write(', normal, ')
				act = 1
			else:
				self.slogf.write(', normal, ')
				act = 0
		else:
			act = self.random_act(rand, Srank)
		self.slogf.write(str(act)+', ')
		self.slogf.close()
		return act


	def select_act_P(self, Q, Srank, Prank):
		rand = random.random()
		if( rand > (self.eps+self.normal) ):
			if( Q[1][Prank] > Q[0][Prank] ):
				self.slogf.write(', normal, ')
				act = 1
			else:
				self.slogf.write(', normal, ')
				act = 0
		else:
			act = self.random_act(rand, Srank)
		self.slogf.write(str(act)+', ')
		self.slogf.close()
		return act


	def select_act_TP(self, Q, Srank, Trank, Prank):
		rand = random.random()
		if( rand > (self.eps+self.normal) ):
			if( Q[1][Trank][Prank] > Q[0][Trank][Prank] ):
				self.slogf.write(', normal, ')
				act = 1
			else:
				self.slogf.write(', normal, ')
				act = 0
		else:
			act = self.random_act(rand, Srank)
		self.slogf.write(str(act)+', ')
		self.slogf.close()
		return act


	def select_act_SP(self, Q, Srank, Prank):
		rand = random.random()
		if( rand > (self.eps+self.normal) ):
			if( Q[1][Srank][Prank] > Q[0][Srank][Prank] ):
				self.slogf.write(', normal, ')
				act = 1
			else:
				self.slogf.write(', normal, ')
				act = 0
		else:
			act = self.random_act(rand, Srank)
		self.slogf.write(str(act)+', ')
		self.slogf.close()
		return act


	def select_act_C(self, Q, Srank, Crank):
		rand = random.random()
		if( rand > (self.eps+self.normal) ):
			if( Q[1][Crank] > Q[0][Crank] ):
				self.slogf.write(', normal, ')
				act = 1
			else:
				self.slogf.write(', normal, ')
				act = 0
		else:
			act = self.random_act(rand, Srank)
		self.slogf.write(str(act)+', ')
		self.slogf.close()
		return act


	def select_act_SC(self, Q, Srank, Crank):
		rand = random.random()
		if( rand > (self.eps+self.normal) ):
			if( Q[1][Srank][Crank] > Q[0][Srank][Crank] ):
				self.slogf.write(', normal, ')
				act = 1
			else:
				self.slogf.write(', normal, ')
				act = 0
		else:
			act = self.random_act(rand, Srank)
		self.slogf.write(str(act)+', ')
		self.slogf.close()
		return act


	def select_act_PC(self, Q, Srank, Prank, Crank):
		rand = random.random()
		if( rand > (self.eps+self.normal) ):
			if( Q[1][Prank][Crank] > Q[0][Prank][Crank] ):
				self.slogf.write(', normal, ')
				act = 1
			else:
				self.slogf.write(', normal, ')
				act = 0
		else:
			act = self.random_act(rand, Srank)
		self.slogf.write(str(act)+', ')
		self.slogf.close()
		return act


	def select_act_SPC(self, Q, Srank, Prank, Crank):
		rand = random.random()
		if( rand > (self.eps+self.normal) ):
			if( Q[1][Srank][Prank][Crank] > Q[0][Srank][Prank][Crank] ):
				self.slogf.write(', normal, ')
				act = 1
			else:
				self.slogf.write(', normal, ')
				act = 0
		else:
			act = self.random_act(rand, Srank)
		self.slogf.write(str(act)+', ')
		self.slogf.close()
		return act



if __name__ == '__main__':
	start = time.time()
	print('start time :', time.ctime())

	action = ['north', 'south', 'east', 'west']
	alpha = 0.1
	gamma = 0.9
	eps = 0.2

	# rewards
	r_move = -1
	r_goal = 100

	# initialize Q table
	Q = np.empty(((4, 4, 4)))


	try:
		for num in range(10000):
			print('----- episode {} -----' .format(num+1))
			state = np.zeros(2, dtype=int)
			nextstate = np.zeros(2, dtype=int)

			while(True):

				# decide action
				if( random.random() > eps ):
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
