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


	def random_act(self):
		if(random.random()>0.5):	heater = 1
		else:						heater = 0
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
		self.interval = interval

		# parameters
		self.alpha = alpha
		self.gamma = gamma

		# epsilon-greedy
		self.eps = 0.05
		self.normal = 0.0

		# log
		self.slogf = open('Qlogs/state_log_'+MODE+'.csv', 'a')
		self.qlogf = open('Qlogs/q_log_'+MODE+'.csv', 'a')


	def initializeQ(self):
		Sr_num = 8
		Tr_num = 6
		Pr_num = 6
		Cr_num = 5
		Rr_num = 7

		# snow accumulation Q[action][Srank]
		if( self.MODE == 'S' ):
			Qtable = np.ones((2, Sr_num))

		# snow accumulation and temperature Q[action][Srank][Trank]
		elif( self.MODE == 'ST' ):
			Qtable = np.ones((2, Sr_num, Tr_num))

		# precipitation Q[action][Prank]
		elif( self.MODE == 'P' ):
			Qtable = np.ones((2, Pr_num))

		# temperature and Precipitation Q[action][Trank][Prank]
		elif( self.MODE == 'TP' ):
			Qtable = np.ones((2, Tr_num, Pr_num))

		# snow accumulation and precipitation Q[action][Srank][Prank]
		elif( self.MODE=='SP' ):
			Qtable = np.ones((2, Sr_num, Pr_num))

		# snow accumulation & continuously operating time
		# Q[action][Srank][Crank]
		elif( self.MODE=='SC' ):
			Qtable = np.ones((2, Sr_num, Cr_num))

		# snow accumulation & road temperature Q[action][Srank][Rrank]
		elif( self.MODE=='SR' ):
			Qtable = np.ones((2, Sr_num, Rr_num))

		# precipitation & road temperature Q[Prank][Rrank]
		elif( self.MODE=='PR' ):
			Qtable = np.ones((2, Pr_num, Rr_num))

		# precipitation and continuously operating time Q[action][Prank][Crank]
		elif( self.MODE=='PC' ):
			Qtable = np.ones((2, Pr_num, Cr_num))

		# snow accumulate, precipitation and continuously operating time
		# Q[action][Srank][Prank][Crank]
		elif( self.MODE=='SPC' ):
			Qtable = np.ones((2, Sr_num, Pr_num, Cr_num))

		# snow accumulate, precipitation and road temperature
		# Q[action][Srank][Prank][Rrank]
		elif( self.MODE=='SPR' ):
			Qtable = np.ones((2, Sr_num, Pr_num, Rr_num))

		elif( self.MODE=='STPR' ):
			Qtable = np.ones((2, Sr_num, Tr_num, Pr_num, Rr_num))

		# invalid MODE
		else:
			print('invalid MODE error')
			sys.exit()

		self.qlogf.write(str(Qtable))
		self.qlogf.close()

		return Qtable


	def Trank(self, temp):
		if( temp < -4 ):	Trank = 0
		elif( temp < -2 ):	Trank = 1
		elif( temp < 0 ):	Trank = 2
		elif( temp < 2 ):	Trank = 3
		elif( temp < 4 ):	Trank = 4
		else:				Trank = 5
		return Trank


	def Srank(self, snow):
		threshold = 0.2
		if( snow < 0 ):
			print('invalid value of snow accumulation')
			sys.exit()
		elif( snow < threshold/2.0 ):	Srank = 0
		elif( snow < threshold ):		Srank = 1
		elif( snow < threshold+0.1 ):	Srank = 2
		elif( snow < threshold+0.2 ):	Srank = 3
		elif( snow < threshold+0.3 ):	Srank = 4
		elif( snow < threshold+1.0 ):	Srank = 5
		elif( snow < 3.0 ):				Srank = 6
		else:							Srank = 7
		return Srank


	def Prank(self, pre):
		if( pre < 0 ):
			print('invalid value of precipitation')
			sys.exit()
		elif( pre == 0 ):	Prank = 0
		elif( pre<0.005 ):	Prank = 1
		elif( pre<0.01 ):	Prank = 2
		elif( pre<0.015 ):	Prank = 3
		elif( pre<0.03 ):	Prank = 4
		else:				Prank = 5
		return Prank

	
	def Crank(self, onT):
		if( onT == 0 ):		Crank = 0
		elif( onT < 10 ):	Crank = 1
		elif( onT < 20 ):	Crank = 2
		elif( onT < 30 ):	Crank = 3
		else:				Crank = 4
		return Crank


	def Rrank(self, TS):
		if( TS < 10 ):		Rrank = 0
		elif( TS < 20 ):	Rrank = 1
		elif( TS < 30 ):	Rrank = 2
		elif( TS < 40 ):	Rrank = 3
		elif( TS < 50 ):	Rrank = 4
		elif( TS < 60 ):	Rrank = 5
		else:				Rrank = 6
		return Rrank


	# calculate rewards
	def calc_reward(self, comp, heater, nextSr):
		reward = 0

		r_on   = -1 * self.interval / 10.0
		r_comp = 100
		r_much = -1 * self.interval / 50.0
		r_nos  = self.interval / 5000.0

		if( nextSr <= 1 ):
			if(comp==False):	reward += r_comp
			if(heater==0):		reward += r_nos
			comp = True
		else:
			if(nextSr>=6 and heater==0):
				reward += r_much
			comp = False

		if(heater==1 and nextSr<7):
			reward += r_on

		return reward, comp


	def nextMax_S(self, Q, eachSr):
		if( Q[1][eachSr] > Q[0][eachSr] ):
			return Q[1][eachSr]
		else:
			return Q[0][eachSr]


	def updateQ_S(self, Q, comp, heater, Srank, onSr, offSr):
		if(heater==1):
			reward, comp = self.calc_reward(comp, heater, onSr)
			nextMax = self.nextMax_S(Q, onSr)
		else:
			reward, comp = self.calc_reward(comp, heater, offSr)
			nextMax = self.nextMax_S(Q, offSr)
		Q[heater][Srank] = (1-self.alpha) * Q[heater][Srank] \
							+ self.alpha * (reward + self.gamma*nextMax)
		Q[heater][Srank] = "%.3f" % Q[heater][Srank]

		self.slogf.write('{:>6}, ' .format(str(reward)))
		self.slogf.close()
		return Q, comp


	def nextMax_ST(self, Q, eachSr, nextTr):
		if( Q[1][eachSr][nextTr] > Q[0][eachSr][nextTr] ):
			return Q[1][eachSr][nextTr]
		else:
			return Q[0][eachSr][nextTr]


	def updateQ_ST(self, Q, comp, heater, Srank, Trank, onSr, offSr, nextTr):
		if(heater==1):
			reward, comp = self.calc_reward(comp, heater, onSr)
			nextMax = self.nextMax_ST(Q, onSr, nextTr)
		else:
			reward, comp = self.calc_reward(comp, heater, offSr)
			nextMax = self.nextMax_ST(Q, offSr, nextTr)
		Q[heater][Srank][Trank] = (1-self.alpha) * Q[heater][Srank][Trank] \
								+ self.alpha * (reward + self.gamma*nextMax)
		Q[heater][Srank][Trank] = '{:.5f}' .format(Q[heater][Srank][Trank])

		self.slogf.write('{:>6}, ' .format(str(reward)))
		self.slogf.close()
		return Q, comp


	def nextMax_P(self, Q, nextPr):
		if( Q[1][nextPr] > Q[0][nextPr] ):
			return Q[1][nextPr]
		else:
			return Q[0][nextPr]


	def updateQ_P(self, Q, comp, heater, onSr, offSr, Prank, nextPr):
		nextMax = self.nextMax_P(Q, nextPr)
		if(heater==1):
			reward, comp = self.calc_reward(comp, heater, onSr)
		else:
			reward, comp = self.calc_reward(comp, heater, offSr)
		Q[heater][Prank] = (1-self.alpha) * Q[heater][Prank]\
							+ self.alpha * (reward + self.gamma*nextMax)
		Q[heater][Prank] = '{:.5f}' .format(Q[heater][Prank])

		self.slogf.write('{:>6}, ' .format(str(reward)))
		self.slogf.close()
		return Q, comp


	def nextMax_TP(self, Q, nextTr, nextPr):
		if( Q[1][nextTr][nextPr] > Q[0][nextTr][nextPr] ):
			return Q[1][nextTr][nextPr]
		else:
			return Q[0][nextTr][nextPr]


	def updateQ_TP(self,Q,comp,heater,onSr,offSr,Trank,Prank,nextTr,nextPr):
		nextMax = self.nextMax_TP(Q, nextTr, nextPr)
		if(heater==1):
			reward, comp = self.calc_reward(comp, heater, onSr)
		else:
			reward, comp = self.calc_reward(comp, heater, offSr)
		Q[heater][Trank][Prank] = (1-self.alpha) * Q[heater][Trank][Prank]\
								+ self.alpha * (reward + self.gamma*nextMax)

		self.slogf.write('{:>6}, ' .format(str(reward)))
		self.slogf.close()
		return Q, comp


	def nextMax_SP(self, Q, eachSr, nextPr):
		if( Q[1][eachSr][nextPr] > Q[0][eachSr][nextPr] ):
			return Q[1][eachSr][nextPr]
		else:
			return Q[0][eachSr][nextPr]


	def updateQ_SP(self, Q, comp, heater, Srank, Prank, onSr, offSr, nextPr):
		if(heater==1):
			reward, comp = self.calc_reward(comp, heater, onSr)
			nextMax = self.nextMax_SP(Q, onSr, nextPr)
		else:
			reward, comp = self.calc_reward(comp, heater, offSr)
			nextMax = self.nextMax_SP(Q, offSr, nextPr)
		Q[heater][Srank][Prank] = (1-self.alpha) * Q[heater][Srank][Prank] \
								+ self.alpha * (reward + self.gamma*nextMax)
		Q[heater][Srank][Prank] = '{:.5f}' .format(Q[heater][Srank][Prank])

		self.slogf.write('{:>6}, ' .format(str(reward)))
		self.slogf.close()
		return Q, comp


	def nextMax_SC(self, Q, eachSr, eachCr):
		if( Q[1][eachSr][eachCr] > Q[0][eachSr][eachCr] ):
			return Q[1][eachSr][eachCr]
		else:
			return Q[0][eachSr][eachCr]


	def updateQ_SC(self, Q, comp, heater, Srank, Crank,\
					onSr, offSr, onCr, offCr):
		if(heater==1):
			reward, comp = self.calc_reward(comp, heater, onSr)
			nextMax = self.nextMax_SC(Q, onSr, onCr)
		else:
			reward, comp = self.calc_reward(comp, heater, offSr)
			nextMax = self.nextMax_SC(Q, offSr, offCr)
		Q[heater][Srank][Crank] = (1-self.alpha) * Q[heater][Srank][Crank]\
									+ self.alpha * (reward + self.gamma*nextMax)
		Q[heater][Srank][Crank] = '{:.5f}' .format(Q[heater][Srank][Crank])

		self.slogf.write('{:>6}, ' .format(str(reward)))
		self.slogf.close()
		return Q, comp


	def nextMax_SR(self, Q, eachSr, eachRr):
		if( Q[1][eachSr][eachRr] > Q[0][eachSr][eachRr] ):
			return Q[1][eachSr][eachRr]
		else:
			return Q[0][eachSr][eachRr]


	def updateQ_SR(self, Q, comp, heater, Srank, Rrank, onSr, offSr, onRr, offRr):
		if(heater==1):
			reward, comp = self.calc_reward(comp, heater, onSr)
			nextMax = self.nextMax_SR(Q, onSr, onRr)
		else:
			reward, comp = self.calc_reward(comp, heater, offSr)
			nextMax = self.nextMax_SR(Q, offSr, offRr)
		Q[heater][Srank][Rrank] = (1-self.alpha)*Q[heater][Srank][Rrank]\
									+ self.alpha*(reward + self.gamma*nextMax)
		Q[heater][Srank][Rrank] = '{:.5f}' .format(Q[heater][Srank][Rrank])

		self.slogf.write('{:>6}, ' .format(str(reward)))
		self.slogf.close()
		return Q, comp


	def nextMax_PR(self, Q, nextPr, eachRr):
		if( Q[1][nextPr][eachRr] > Q[0][nextPr][eachRr] ):
			return Q[1][nextPr][eachRr]
		else:
			return Q[0][nextPr][eachRr]


	def updateQ_PR(self,Q,comp,heater,Prank,Rrank,onSr,offSr,nextPr,onRr,offRr):
		if(heater==1):
			reward, comp = self.calc_reward(comp, heater, onSr)
			nextMax = self.nextMax_PR(Q, nextPr, onRr)
		else:
			reward, comp = self.calc_reward(comp, heater, offSr)
			nextMax = self.nextMax_PR(Q, nextPr, offRr)
		Q[heater][Prank][Rrank] = (1-self.alpha)*Q[heater][Prank][Rrank]\
									+ self.alpha*(reward + self.gamma*nextMax)
		Q[heater][Prank][Rrank] = '{:.5f}' .format(Q[heater][Prank][Rrank])

		self.slogf.write('{:>6}, ' .format(str(reward)))
		self.slogf.close()
		return Q, comp


	def nextMax_PC(self, Q, nextPr, nextCr):
		if( Q[1][nextPr][nextCr] > Q[0][nextPr][nextCr] ):
			return Q[1][nextPr][nextCr]
		else:
			return Q[0][nextPr][nextCr]


	def updateQ_PC(self,Q,comp,heater,onSr,offSr,Prank,Crank,nextPr,onCr,offCr):
		if(heater==1):
			reward, comp = self.calc_reward(comp, heater, onSr)
			nextMax = self.nextMax_PC(Q, nextPr, onCr)
		else:
			reward, comp = self.calc_reward(comp, heater, offSr)
			nextMax = self.nextMax_PC(Q, nextPr, offCr)
		Q[heater][Prank][Crank] = (1-self.alpha) * Q[heater][Prank][Crank]\
								+ self.alpha * (reward + self.gamma*nextMax)
		Q[heater][Prank][Crank] = '{:.5f}' .format(Q[heater][Prank][Crank])

		self.slogf.write('{:>6}, ' .format(str(reward)))
		self.slogf.close()
		return Q, comp


	def nextMax_SPC(self, Q, nextSr, nextPr, nextCr):
		if( Q[1][nextSr][nextPr][nextCr] > Q[0][nextSr][nextPr][nextCr] ):
			return Q[1][nextSr][nextPr][nextCr]
		else:
			return Q[0][nextSr][nextPr][nextCr]


	def updateQ_SPC(self, Q, comp, heater, Srank, Prank, Crank, onSr, offSr,\
					nextPr, onCr, offCr):
		if(heater==1):
			reward, comp = self.calc_reward(comp, heater, onSr)
			nextMax = self.nextMax_SPC(Q, onSr, nextPr, onCr)
		else:
			reward, comp = self.calc_reward(comp, heater, offSr)
			nextMax = self.nextMax_SPC(Q, offSr, nextPr, offCr)
		h = heater
		Q[h][Srank][Prank][Crank] = (1-self.alpha) * Q[h][Srank][Prank][Crank]\
									+ self.alpha * (reward + self.gamma*nextMax)
		Q[h][Srank][Prank][Crank] = '{:.2f}' .format(Q[h][Srank][Prank][Crank])

		self.slogf.write('{:>6}, ' .format(str(reward)))
		self.slogf.close()
		return Q, comp


	def nextMax_SPR(self, Q, nextSr, nextPr, nextRr):
		if( Q[1][nextSr][nextPr][nextRr] > Q[0][nextSr][nextPr][nextRr] ):
			return Q[1][nextSr][nextPr][nextRr]
		else:
			return Q[0][nextSr][nextPr][nextRr]


	def updateQ_SPR(self, Q, comp, heater, Srank, Prank, Rrank, onSr, offSr,\
					nextPr, onRr, offRr):
		if(heater==1):
			reward, comp = self.calc_reward(comp, heater, onSr)
			nextMax = self.nextMax_SPR(Q, onSr, nextPr, onRr)
		else:
			reward, comp = self.calc_reward(comp, heater, offSr)
			nextMax = self.nextMax_SPR(Q, offSr, nextPr, offRr)
		h = heater
		Q[h][Srank][Prank][Rrank] = (1-self.alpha) * Q[h][Srank][Prank][Rrank]\
									+ self.alpha * (reward + self.gamma*nextMax)
		Q[h][Srank][Prank][Rrank] = '{:.2f}' .format(Q[h][Srank][Prank][Rrank])

		self.slogf.write('{:>6}, ' .format(str(reward)))
		self.slogf.close()
		return Q, comp


	def nextMax_STPR(self, Q, nextSr, nextTr, nextPr, nextRr):
		if( Q[1][nextSr][nextTr][nextPr][nextRr]\
					> Q[0][nextSr][nextTr][nextPr][nextRr] ):
			return Q[1][nextSr][nextTr][nextPr][nextRr]
		else:
			return Q[0][nextSr][nextTr][nextPr][nextRr]


	def updateQ_STPR(self, Q, comp, heater, Srank, Trank, Prank, Rrank, onSr,\
							offSr, nextTr, nextPr, onRr, offRr):
		if(heater==1):
			reward, comp = self.calc_reward(comp, heater, onSr)
			nextMax = self.nextMax_STPR(Q, onSr, nextTr, nextPr, onRr)
		else:
			reward, comp = self.calc_reward(comp, heater, offSr)
			nextMax = self.nextMax_STPR(Q, offSr, nextTr, nextPr, offRr)
		h = heater
		Q[h][Srank][Trank][Prank][Rrank] = (1-self.alpha)\
											* Q[h][Srank][Trank][Prank][Rrank]\
										+ self.alpha*(reward + self.gamma*nextMax)
		Q[h][Srank][Trank][Prank][Rrank] = '{:.2f}'\
										.format(Q[h][Srank][Trank][Prank][Rrank])

		self.slogf.write('{:>6}, ' .format(str(reward)))
		self.slogf.close()
		return Q, comp


	def random_act(self, rand, Srank):
		if( rand <= (self.eps+self.normal) and rand > self.eps ):
			self.slogf.write(', n, ')
			if(Srank==0):		act = 0
			else:				act = 1
		else:
			self.slogf.write(', r, ')
			act = random.choice(action)
		return act


	def select_act_S(self, Q, Srank):
		rand = random.random()
		if( rand > (self.eps+self.normal) ):
			if( Q[1][Srank] >= Q[0][Srank] ):
				self.slogf.write(', n, ')
				act = 1
			else:
				self.slogf.write(', n, ')
				act = 0
		else:
			act = self.random_act(rand, Srank)
		self.slogf.write(str(act)+', ')
		self.slogf.close()
		return act


	def select_act_ST(self, Q, Srank, Trank):
		rand = random.random()
		if( rand > (self.eps+self.normal) ):
			if( Q[1][Srank][Trank] >= Q[0][Srank][Trank] ):
				self.slogf.write(', n, ')
				act = 1
			else:
				self.slogf.write(', n, ')
				act = 0
		else:
			act = self.random_act(rand, Srank)
		self.slogf.write(str(act)+', ')
		self.slogf.close()
		return act


	def select_act_P(self, Q, Srank, Prank):
		rand = random.random()
		if( rand > (self.eps+self.normal) ):
			if( Q[1][Prank] >= Q[0][Prank] ):
				self.slogf.write(', n, ')
				act = 1
			else:
				self.slogf.write(', n, ')
				act = 0
		else:
			act = self.random_act(rand, Srank)
		self.slogf.write(str(act)+', ')
		self.slogf.close()
		return act


	def select_act_TP(self, Q, Srank, Trank, Prank):
		rand = random.random()
		if( rand > (self.eps+self.normal) ):
			if( Q[1][Trank][Prank] >= Q[0][Trank][Prank] ):
				self.slogf.write(', n, ')
				act = 1
			else:
				self.slogf.write(', n, ')
				act = 0
		else:
			act = self.random_act(rand, Srank)
		self.slogf.write(str(act)+', ')
		self.slogf.close()
		return act


	def select_act_SP(self, Q, Srank, Prank):
		rand = random.random()
		if( rand > (self.eps+self.normal) ):
			if( Q[1][Srank][Prank] >= Q[0][Srank][Prank] ):
				self.slogf.write(', n, ')
				act = 1
			else:
				self.slogf.write(', n, ')
				act = 0
		else:
			act = self.random_act(rand, Srank)
		self.slogf.write(str(act)+', ')
		self.slogf.close()
		return act


	def select_act_SC(self, Q, Srank, Crank):
		rand = random.random()
		if( rand > (self.eps+self.normal) ):
			if( Q[1][Srank][Crank] >= Q[0][Srank][Crank] ):
				self.slogf.write(', n, ')
				act = 1
			else:
				self.slogf.write(', n, ')
				act = 0
		else:
			act = self.random_act(rand, Srank)
		self.slogf.write(str(act)+', ')
		self.slogf.close()
		return act


	def select_act_SR(self, Q, Srank, Rrank):
		rand = random.random()
		if( rand > (self.eps+self.normal) ):
			if( Q[1][Srank][Rrank] >= Q[0][Srank][Rrank] ):
				self.slogf.write(', n, ')
				act = 1
			else:
				self.slogf.write(', n, ')
				act = 0
		else:
			act = self.random_act(rand, Srank)
		self.slogf.write(str(act)+', ')
		self.slogf.close()
		return act


	def select_act_PR(self, Q, Srank, Prank, Rrank):
		rand = random.random()
		if( rand > (self.eps+self.normal) ):
			if( Q[1][Prank][Rrank] >= Q[0][Prank][Rrank] ):
				self.slogf.write(', n, ')
				act = 1
			else:
				self.slogf.write(', n, ')
				act = 0
		else:
			act = self.random_act(rand, Srank)
		self.slogf.write(str(act)+', ')
		self.slogf.close()
		return act

	def select_act_PC(self, Q, Srank, Prank, Crank):
		rand = random.random()
		if( rand > (self.eps+self.normal) ):
			if( Q[1][Prank][Crank] >= Q[0][Prank][Crank] ):
				self.slogf.write(', n, ')
				act = 1
			else:
				self.slogf.write(', n, ')
				act = 0
		else:
			act = self.random_act(rand, Srank)
		self.slogf.write(str(act)+', ')
		self.slogf.close()
		return act


	def select_act_SPC(self, Q, Srank, Prank, Crank):
		rand = random.random()
		if( rand > (self.eps+self.normal) ):
			if( Q[1][Srank][Prank][Crank] >= Q[0][Srank][Prank][Crank] ):
				self.slogf.write(', n, ')
				act = 1
			else:
				self.slogf.write(', n, ')
				act = 0
		else:
			act = self.random_act(rand, Srank)
		self.slogf.write(str(act)+', ')
		self.slogf.close()
		return act


	def select_act_SPR(self, Q, Srank, Prank, Rrank):
		rand = random.random()
		if( rand > (self.eps+self.normal) ):
			if( Q[1][Srank][Prank][Rrank] >= Q[0][Srank][Prank][Rrank] ):
				self.slogf.write(', n, ')
				act = 1
			else:
				self.slogf.write(', n, ')
				act = 0
		else:
			act = self.random_act(rand, Srank)
		self.slogf.write(str(act)+', ')
		self.slogf.close()
		return act


	def select_act_STPR(self, Q, Srank, Trank, Prank, Rrank):
		rand = random.random()
		if( rand > (self.eps+self.normal) ):
			if(Q[1][Srank][Trank][Prank][Rrank]>=Q[0][Srank][Trank][Prank][Rrank]):
				self.slogf.write(', n, ')
				act = 1
			else:
				self.slogf.write(', n, ')
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
	Q = np.ones(((4, 4, 4)))


	try:
		for num in range(10000):
			print('----- episode {} -----' .format(num+1))
			state = np.ones(2, dtype=int)
			nextstate = np.ones(2, dtype=int)

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
