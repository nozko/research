#!/usr/bin/env python
#-*- coding: utf-8 -*-

from __future__ import print_function
import math
import numpy as np
import re
import math

RD = 3.1416/180
pido = 43.04
pked = 141.33
pcth = 0
pctl = 0.2
tvm = (pcth + pctl)/2
cth = 1.67


class sim:

	def __init__(self):
		fn1 = 'file.net'
		self.net = open(fn1, 'r')
		fn2 = 'op.csv'
		self.op = open(fn2, 'a')
		fn3 = 'smap.run'
		self.run = open(fn3, 'r')
		fn4 = 'sekisetu.csv'
		self.snow = open(fn4, 'a')
		fn5 = 'asahikawa.csv'
		self.data = open(fn5, 'r')
		fn6 = 'hist.rst'
		self.hist = open(fn6, 'r')
		fn7 = 'file2.net'
		self.net2 = open(fn7, 'r')


	# absolute humidity
	def abshumid(self, Temp):
		P = 51.7063 * 10**(6.21147-2886.37/(1.8*Temp+491.69)\
				- 337269.46/(1.8*Temp+491.69)**2)
		funx = 0.622*P / (760-P)
		return funx


	def funa(self, Wspeed):
		funa = 5 + 3.4*Wspeed
		return funa

	
	def sunpo(self, TM, SH, CH, SA, CA):
		W = 2 * 3.1416 * Lday / 366
		Del = 0.362213 - 23.2476*math.cos(W+0.153231) \
			- 0.336891*math.cos(2*W+0.207099) - 0.185265*math.cos(3*W+0.620129)
		AE = -0.000278641 + 0.122772*math.cos(W+1.49831) \
			- 0.165458*math.cos(2*W-1.26155) - 0.00535383*math.cos(3*W-1.1571)
		sjik = 15 * (TM-12+AE) + pked - 135
		SH = math.sin(pido*RD)*math.sin(Del*RH) \
			+ math.cos(pido*RD)*math.cos(Del*RD)*math.cos(sjik*RD)
		if(SH<=0):
			SH = 0
			CH = 0
			SA = 0
			CA = 0
		else:
			CH = math.sqrt(1-SH**2)
			SA = math.cos(Del*RD)*math.sin(sjik*RD)/CH
			CA = (SH*math.sin(pido*RD)-math.sin(Del*RD)) / CH / math.cos(pido*RD)
		return SH, CH, SA, CA


	def qds(SH, CH, SA, CA, rdn, rsh, rtv):
		if(SH > 0):
			rdv = rdn * SH
		else:
			rdv = 0
		rfv = rsh
		rtv = rdv + rfv
		return rdv, rfv, rtv



if __name__ == '__main__':

	# initialize array
	BT = []
	Q = []
	E = []
	tset = []
	ID = []
	QQ = []
	npn = []
	nsd = []
	nqd = []
	C = []
	CR = []
	AS = []
	iopt = []
	iope = []
	iopq = []
	CL = []
	ict = []
	TD = []
	QD = np.zeros(100)
	ED = np.zeros(100)
	lms = []
	lds = []
	jdn = []
	jtd = []
	jqd = []
	jed = []
	QE = []
	BF = []
	htr = []
	NP = []
	nkr = []
	HR = []
	S1 = []
	S2 = []
	TT = []
	ER = []
	flr = []
	IL = []
	BW = np.zeros(5000)
	BS = []
	SS = []
	WW = []
	EV = []
	Scover = []
	bdph = []
	mlt = []
	sat = []
	NC = []
	smf = []
	Theight = []
	Tcover = []
	nfd = np.zeros(100)
	nsca = np.zeros(100)
	dpd = np.empty(100)
	ssd = np.empty(100)
	tsd = np.empty(100)


	# initialize
	Lyear  = 0
	temp_o = 0		# temperature
	Qsup   = 0		# supply calorie
	RT     = 0		# operating time
	height = 0		# height of snow
	slevel = 0		# level of snow accumulation
	etsx   = 0
	idr    = 2
	irt    = 0
	melt   = 0		# melted snow
	plt    = 0
	tmr    = 0
	lpmx   = 0
	idx    = int( 1/(interval+0.4) )
	lpr    = 1.0
	pct    = 0
	ntime  = np.zeros((5, 5))
	Qr  = 0		# Q from surface(snowing, plus temperature)
	Qr2 = 0		# Q from surface(snowing, plus temperature) @moist sensor node


	### read file 'smap.run'
	smap1 = sim().run.readline().split('\t')
	interval   = float(smap1[0])	# calculatioin interval
	outputStep = float(smap1[1])	# output step interval
	diffM      = int(smap1[2])		# difference method(1:後退差分)
	CK         = float(smap1[3])	# 緩和係数(0.7~1.5程度)
	maxloop    = int(smap1[4])		# maximum numper of loops(200)
	print('interval :', interval, outputStep, diffM, CK, maxloop)

	smap2 = sim().run.readlines()[1].split('\t')
	Dstart = int(smap2[0])	# day to start calculation(days from top of the file)
	lyear  = int(smap2[1])	# number of years to calculate(1~)
	Dend   = int(smap2[2])	# last day of calculation(days from top of the file)
	print('start :', Dstart, 'years :', lyear, 'end :', Dend)
	Lday = Dstart - 1

	Dop = int( sim().run.readlines()[2] )	# days to output in detail('op.csv')
	for i in range(Dop):
		smaps = sim().run.readlines()[3+i].split('\t')
		lms.append(smaps[0])
		lds.append(smaps[1])

	smap4 = sim().run.readlines()[4+Dop].split('\t')
	Qs      = float(smap4[0])	# Q of heat source
	flowR   = float(smap4[1])	# flow rate
	maxT    = float(smap4[2])	# maximum temperature of heat source
	Tsensor = int(smap4[3])		# number of node(temperature sensor)
	Wsensor = int(smap4[4])		# number of node(moist sensor)
	timer   = float(smap4[5])	# delay time
	step    = int(smap4[6])		# rotation step
	circuit = int(smap4[7])		# number of rotatioin circuit
	print(Qs, flowR, maxT, Tsensor, Wsensor, timer, step, circuit)

	smap5 = sim().run.readlines()[5+Dop].split('\t')
	snowT  = float(smap5[0])	# operating temperature during snowfall
	B1     = float(smap5[1])	# 外気温にかかる係数1
	stopT1 = float(smap5[2])	# stop temperature(difference from operating one)
	wetT   = float(smap5[3])	# operating temperature with wet surface
	B2     = float(smap5[4])	# 外気温にかかる係数2
	stopT2 = float(smap5[5])	# stop temperature
	dryT   = float(smap5[6])	# operating temperature with dry surface
	B3     = float(smap5[7])	# 外気温にかかる係数3
	stopT3 = float(smap5[8])	# stop temperature
	offT   = float(smap5[9])	# not driving temperature
	print(snowT, B1, stopT1, wetT, B2, stopT2, dryT, B3, stopT3, offT)

	smap6 = sim().run.readlines()[6+Dop].split('\t')
	remainW = float(smap4[0])	# remained water after melting
	maxPene = float(smap4[1])	# maximum penetration height of water
	abrate0 = float(smap4[2])	# solar radiatioin absorption rate of snow(0.2)
	windC   = float(smap4[3])	# wind speed correction coefficient(~1)
	Scover0 = float(smap4[4])	# initial snow cover
	dens0   = float(smap4[5])	# initial snow density
	BS0 = Scover0 * dens0

	sim().run.close()


	# count initialize
	itn = 0
	iqn = 0
	ien = 0
	idn = 0


	line_nun = 0
	ipx = int(sim().net.readlines()[line_num])		# data num

	### read all data in 'file.net'
	for ip in range(ipx):
		line_num += 1
		net1 = re.split(" +", sim().net.readlines()[line_num])
		print(net1)
		DM = net1[1]
		BT.append(float(net1[2]))
		ID.append(float(net1[3]))
		C.append(float(net1[4]))
		CL.append(float(net1[5]))

		line_num += 1
		net2 = re.split(" +", sim().net.readlines()[line_num])
		print(net2)
		nsd.append(float(net2[1]))
		AS.append(float(net2[2]))
		nqd.append(float(net2[3]))
		QQ.append(float(net2[4]))

		Scover.append(0)
		BS.append(0)

		line_num += 1
		net3 = re.split(" +", sim().net.readlines()[line_num])
		print(net3)
		iopt.append(int(net3[1]))
		iopq.append(float(net3[2]))
		iope.append(float(net3[3]))
		if(int(net3[1])==1):
			itn += 1
			jtd.append(ip)
			iqn += 1
			jqd.append(ip)
			ien += 1
			jed.append(ip)
		if(int(net2[1])==1):
			area = area + float(net2[2])
			BS.append(BS0)
			Scover.append(Scover0)
			bdph.append(Scover0)
			idn += 1
			jdn.append(ip)

		line_num += 1
		net4 = re.split(" +", sim().net.readlines()[line_num])
		print(net4)
		npn.append(int(net4[1]))
		_NP = []
		_nkr = []
		_HR = []
		for j in range(int(net4[1])):
			line_num += 1
			net4j = re.split(" +", sim().net.readlines()[line_num])
			_NP.append(int(net4j[1]))
			_nkr.append(float(net4j[2]))
			_HR.append(float(net4j[3]))
		NP.append(_NP)
		nkr.append(_nkr)
		HR.append(_HR)

		flr.append(0)

	sim().net.close()


	ih = ipx + 1


	# output 'op.csv'
	sim().op.write('month,day,TM,temp_o,absH,sun,nightR,Wspeed,pre,\
						FD,sca,water,evaporate,SM,SE,flow')
	for i in range(itn):
		sim().op.write(', ' + str(jtd[i]))
	for i in range(iqn):
		sim().op.write(', ' + str(jqd[i]))
	for i in range(ien):
		sim().op.write(', ' + str(jed[i]))
	sim().op.write('\n')
	
	sim().snow.write('month,day,TM')
	for i in range(idn):
		for j in range(3):
			sim().snow.write(', ' + str(jdn[1]))
	sim().snow.write('\n')
	

	### day loop ###
	for day in range(Dend-Dstart):
		Lday += 1
		if(Lday > 365):	Lday = 1
		print(Lday, Dstart)
		if(Lday==Dstart):	Lyear += 1
		print(Lyear, Lday)


		ids = 0
		for i in range(Dop):
			if( MT==lms[i] and MD==lds[i] ):
				ids = 1


		### time loop ###
		for t in range(24):
			t0 = temp_o
			print(temp_o)
			data1 = sim().data.readlines()[t].split(',')
			month   = int(data1[1])
			day     = int(data1[2])
			tempA   = float(data1[4])	# temperature
			VP      = float(data1[5])
			Wspeed0 = float(data1[6])
			sun     = float(data1[7])	# solar radiatioin
			pre     = float(data1[8])	# precipitation降水量
			unr     = float(data1[9])
			nightR  = 45				# nighttime radiation

			sun = sun/4.186*1000

			if(VP >= 0):
				absH = 0.622*VP/(1013.25-VP)	# absolute humidity
				RH = VP*760/1013.25
				if(unr >= 0):
					nightR = 0.0000000488 * (tempA+273.16)**4 \
							* (1-0.62*unr/10) * (0.49-0.076*math.sqrt(RH))
			else:
				absH = 0.7 * sim().abshumid(tempA)

			for i in range(Dop):
				if( month==lms[i] and day==lds[i] ):
					ids = 1

			for idt in range(idx):
				TM = t - 1 + interval*idt
	
				temp_o = t0 + (tempA-t0)*interval*idt
				Wspeed = Wspeed0 * windC
				if(temp_o < 0):
					FS = pre
					FR = 0
				elif(temp_o>=0 and temp_o<2):
					FS = pre * (2-temp_o)/2
					FR = pre * temp_o / 2
				elif(temp_o>=2):
					FS = 0
					FR = pre

				YG = 1000 / (0.091*temp_o**2 - 1.81*temp_o + 9.47)
				if(YG < 50):	YG = 50


				# set
				sca = 0
				ssmf = 0
				FD = 0
				water = 0		# water amount
				SE = 0			# heat source calorific value
				erot = 0
				evaporate = 0	# evaporation amount
				ilp = 0
				erx = 0
				L1 = 0
				L2 = 0
				if(idr==2):
					pmx = 0
					ptl = 0
					ptm = 0
				if(pre > 0):
					ptl = plt + pre*interval
					ptm = ptm + interval
				if(pre > pmx):
					pmx = pre
				if(snowT==99):
					wetT = pmx * 2
					if(wetT < 2):	wetT = 2
					B1 = wetT
					B2 = wetT
				if(snowT==999):
					wetT = pmx + 2
					dryT = wetT
					B1 = wetT
					B2 = wetT
					B3 = wetT
				tmr = tmr + interval

				if(pre > 0):	tmr = 0
				if(month<=4 or month>=11):
					if(pre > 0):
						level = 1
					else:
						if( (BS[Wsensor]+BW[Wsensor]) > 0 ):
							level = 2
						else:
							level = 3
					if(level==1):
						TN = snowT + B1*temp_o
						TF = TN + stopT1
					elif(level==2):
						TN = wetT + B2*temp_o
						TF = TN + stopT2
					elif(level==3):
						TN = S3 + B3*temp_o
						TF = TN + stopT3

					if(BT[Tsensor] < TN):		idr = 1
					elif(BT[Tsensor] > TF):		idr = 2
					if(tmr > timer+0.0001):		idr = 2
					if(temp_o > offT):			idr = 2

					if(idr==2):
						irt = -1
						imo = 0
					elif(idr==1):
						irt += 1
						IA = irt % (circuit*step)	# circuit*step : rotation time
						if(IA>=0 and IA<step):			imo = 1
						elif(IA>=step and IA<step*2):	imo = 2
						elif(IA>=step*2 and IA<step*3):	imo = 3
					ntime[level][idr] += 1

				for ip in range(ipx):
					E.append(0)
					Q.append(0)
					NC.append(0)
					mlt.append(0)
					ict.append(0)

					if(ID[ip]==2):
						tset.append(temp_o)
						ict[ip] = 1
					elif(ID[ip]==9):
						tset.append(BT[ip])
						ict[ip] = 1

					if(nqd[ip]==1):
						Q[ip] = QQ[ip]
					elif(nqd[ip]==11):
						if(idr==1):
							ih = ip
							Q[ip] = Qs
							ict[ip] = 0

					if(idr==1):
						if(ID[ip]==10):		# normal heater
							Q[ip] = Qs
							erot = erot + Q[ip]
						elif(ID[ip]>=11 and ID[ip]<=13):
							Q[ip] = Qs * 0.5
							if(imo>=1 and imo<=3):	Q[ip] = Qs*2
							erot = erot + Q[ip]

					for j in range(npn[ip]):
						if(nkr[ip][j]==11):
							HR[ip][j] = 0
							if(idr==1):
								HR[ip][j] = flowR
							flow = HR[ip][j]

					CR.append(C[ip])
					if(CL[ip] > 0):
						if(BT[ip] > tvm):
							CR[ip] = C[ip] + CL[ip]
						else:
							CR[ip] = C[ip] + CL[ip]*0.5
						Q[ip] = flr[ip]
						if(IL[ip]==1):
							tic = BT[ip]
							pic = (pcth-tic) / (pcth-pctl)
							CR[ip] = C[ip] + CL[ip]*(1-pic) \
									+ CL[ip]*0.5*pic + CL[ip]*80/abs(pcth-pctl)

					if(nsd[ip]==1):
						abrate = 0.8 - 30*bdph[ip]
						if(abrate < abrate0):	abrate = abrate0
						sat.append( temp_o \
								+ (abrate*sun-0.9*nightR)/(sim().funa(Wspeed)+4) )
						EV.append(0)
						TS = BT[ip]
						lps = 0
						dps = bdph[ip]
						if( (BS[ip]+FS)>0 ):
							if(TS<0 or sat[ip]>0):
								mlt[ip] = 1
						if( (BW[ip]+FR)>0 ):
							if(TS < 0):
								mlp[ip] = 1
						if( (BS[ip]+FS)>0 and (BW[ip]+FR)>0 ):
							mlt[ip] = 1
						HI = 0
						if(bdph[ip] > 0):
							HI = BW[ip] / (1000 - BS[ip]/bdph[ip])
						else:
							HI = BW[ip] / 1000
						if(HI > maxPene):
							HI = maxPene
						wat = BW[ip] + FR*interval

# 49
						SM = 0
						EV[ip] = 0
						BF.append(0)
						if(mlt[ip]==1):
							DH = dps = HI
							if(DH<=0 or sat[ip]>0):
								DH = 0
								tsv = 0
								EV[ip] = 4 * sim().funa(Wspeed) \
										* (sim().abshumid(tsv)-absH)
								wat = BW[ip] + (FR-EV[ip])*interval
								if(wat < 0):
									EV[ip] = BW[ip]/interval + FR
									wat = 0
							htrm = 1 / (1/(sim().funa(Wspeed)+4) + DH/0.08)
							SM0 = (200*TS + htrm*sat[ip] - 590*EV[ip])*interval/80
							SM = SM0
							BF[ip] = 1
							if(SM0 < -1*wat):
								SM = -1 * wat
								BF[ip] = SM / SM0
							elif( SM0 > (BS[ip]+FS*interval) ):
								SM = BS[ip] + FS*interval
								BF[ip] = SM / SM0
						elif(mlt[ip]==0):
							tsv = BT[ip]
							EV[ip] = 4 * sim().funa(Wspeed)\
									* (sim().abshumid(tsv)- absH)
							if( EV[ip] > (BW[ip]/interval + FR) ):
								EV[ip] = BW[ip]/interval + FR

						htr.append( 1 / (1/(sim().funa(Wspeed)+4) + dps/0.08) )
						QE.append( -590*EV[ip] * AS[ip] * (1-BF[ip]) )
						SS.append( BS[ip] - SM + FS*interval )
						WW.append( BW[ip] + SM + (FR-EV[ip])*interval )
						smf.append( SM )

						if(SS[ip] > 0):
							if(SM > 0):
								G = (BS[ip] + FS*interval) \
									/ (bdph[ip] + FS*interval/YG)
								Scover[ip] = SS[ip] / G
							else:
								Scover[ip] = bdph[ip] + FS*interval/YG - SM/916
						else:
							Scover[ip] = 0

						ddps = (Scover[ip]+bdph[ip])/2 - dps
						if(abs(ddps)>0.001 and lps<10):
							dps = 0.7*ddps + dps
							lps += 1
							# go to 49


				# matrix
				T = BT
# 71
				lps = 0
# 76
				for ip in range(ipx):
					if(ict[ip]==1):
						T[ip] = tset[ip]
						E[ip] = 0
					S1.append( (Q[ip]+QE[ip]+E[ip])*interval + BT[ip]*CR[ip] )
					S2.append( CR[ip] )

					for j in range(npn[ip]):
						np_n = NP[ip][j]
						tmp = diffM*T[np_n] + (1-diffM)*(BT[np_n]-BT[ip])
						S1[ip] = S1[ip] + HR[ip][j]*(tmp)*interval
						S2[ip] = S2[ip] + diffM*HR[ip][j]*interval

					if(nsd[ip]==1):
						F = BF[ip]
						tmp = (1-F)*htr[ip]*sat[ip] - F*200*(1-diffM)*BT[ip] \
								- (1-F)*htr[ip]*(1-diffM)*BT[ip]
						S1[ip] = S1[ip] + ( tmp )*interval*AS[ip]
						S2[ip] = S2[ip] \
								+ (F*200+(1-F)*htr[ip])*interval*AS[ip]*diffM
					if(ict[ip]==1):
						E[ip] = (S2[ip]*tset[ip] - S1[ip]) / interval
					else:
						if(diffM<=0.01):
							T[ip] = S1[ip] / S2[ip]
							erx = 0
						else:
							TT.append( S1[ip] / S2[ip] )
							ER.append( TT[ip] - T[ip] )
							aer = abs(ER[ip])
							if(aer > erx):
								erx = aer
								ier = ip
							T[ip] = CK*ER[ip] + T[ip]

				lps += 1
				if(lps > maxloop):
					print(month, day, TM, erx, ier, T[ier])
				else:
					if(erx > 0.0001):
						erx = 0
						# go to 76
					if(lps > lpmx):	lpmx = lps

				if(ilp < 2):
					if( T[ih]>maxT and Q[ih]>0 ):
						ilp += 1
						tset[ih] = maxT
						ict[ih] = 1
						Q[ih] = 0
						erx = 0
						# go to 71
					if(E[ih] < 0):
						ilp += 1
						ict[ih] = 0
						E[ih] = 0
						erx = 0
						# go to 71

				for ip in range(ipx):
					if(CL[ip] > 0):
						IL.append(0)
						flr[ip] = 0
						if( T[ip]<pcth and T[ip]>pctl ):
							IL[ip] = 1
						elif( BT[ip]>pcth and T[ip]<pcth ):
							flr[ip] = (C[ip]+CL[ip]) * (T[ip]-pcth)
							T[ip] = pcth
							IL[ip] = 1
						elif( BT[ip]>pctl and T[ip]<pctl ):
							flr[ip] = (C[ip] + CL[ip]*80/abs(pcth-pctl)) \
										* (T[ip]-pctl)
							T[ip] = pctl
						elif( BT[ip]<pctl and T[ip]>pctl ):
							flr[ip] = (C[ip] + CL[ip]*0.5) * (T[ip]-pctl)
							T[ip] = pctl
							IL[ip] = 1
						elif( BT[ip]<pcth and T[ip]>pcth ):
							flr[ip] = (C[ip] + CL[ip]*80/abs(pcth-pctl)) \
										* (T[ip]-pcth)
							T[ip] = pcth

					if(nsd[ip]==1):
						if( level==1 and T[ip]<0 ):
							L1 = 1
						elif( level==2 and T[ip]<0 ):
							L2 = 1
						if(SS[ip] < 0):
							SS[ip] = 0
						if(WW[ip] < 0):
							WW[ip] = 0
						if(Scover[ip] > FD):
							FD = Scover[ip]

						sca = sca + SS[ip]*AS[ip]/area
						water = water + WW[ip]*AS[ip]/area
						evaporate = evaporate + EV[ip]*interval*AS[ip]/area
						ssmf = ssmf + smf[ip]*AS[ip]/area

						if( Scover[ip]>0 and SS[ip]>0 ):
							GM = SS[ip] / Scover[ip]
							EE = 16 * math.exp(0.021*GM)
							GM = GM * math.exp( SS[ip]/2/EE*interval/24 )
							if(GM > 916):
								GM = 916
							Scover[ip] = SS[ip] / GM
				

				if(FD > 0):
					ifd = FD * 100
					for k in range(ifd+1):
						nfd[k] += 1
					if(ifd > height):
						height = ifd
						print(height, ifd, fd)

				if(sca > 0):
					isca = sca
					for k in range(isca+1):
						nsca[k] += 1
					if(isca > slevel):
						slevel = isca

				snow_minusT += L1
				wet_minusT  += L2
				SE = E[ih] + Q[ih] + erot

				if(SE > 0):
					RT += 1
					Qsup += se
					rse = Qsup / RT

				melt += ssmf

				if(ids>=0):
					if(lpr==outputStep):
						# count initialize
						itn = 0
						iqn = 0
						ien = 0
						idp = 0

						for ip in range(ipx):
							if(iopt[ip]==1):
								itn += 1
								TD[itn] = T[ip]
							if(iopq[ip]==1):
								iqn += 1
								QD[iqn] = Q[ip]
							if(iope[ip]==1):
								ien += 1
								ED[ien] = E[ip]
							if(nsd[ip]==1):
								idp += 1
								dpd[idp] = Scover[ip]
								ssd[idp] = SS[ip]
								tsd[idp] = T[ip]

						sim().op.write(str(month) + ', ')
						sim().op.write(str(day) + ', ')
						sim().op.write(str(TM) + ', ')
						sim().op.write(str(temp_o) + ', ')
						sim().op.write(str(absH) + ', ')
						sim().op.write(str(sun/0.86) + ', ')
						sim().op.write(str(nightR/0.86) + ', ')
						sim().op.write(str(Wspeed) + ', ')
						sim().op.write(str(pre) + ', ')
						sim().op.write(str(FD) + ', ')
						sim().op.write(str(sca) + ', ')
						sim().op.write(str(water) + ', ')
						sim().op.write(str(evaporate) + ', ')
						sim().op.write(str(ssmf) + ', ')
						sim().op.write(str(SE/0.86/area) + ', ')
						sim().op.write(str(flow) + ', ')
						for i in range(itn):
							sim().op.write(str(TD[i]) + ', ')
						for i in range(iqn):
							sim().op.write(str(QD[i]/0.86) + ', ')
						for i in range(ien):
							sim().op.write(str(ED[i]/0.86) + ', ')
						sim().op.write('\n')

						sim().snow.write(str(month) + ', ')
						sim().snow.write(str(day) + ', ')
						sim().snow.write(str(TM) + ', ')
						for i in range(idp):
							sim().snow.write(str(dpd[i]) + ', ')
						for i in range(idp):
							sim().snow.write(str(ssd[i]) + ', ')
						for i in range(idp):
							sim().snow.write(str(tsd[i]) + ', ')

						lpr = 0

					lpr += 1

				if( pre>0 and T[Is]>0 ):
					TS = diffM*T[Wsensor] + (1-diffM)*BT[Wsensor]
					Qr2 = Qr2 + BF[Wsensor]*200*TS \
							+ (1-BF[Wsensor])*htr[Wsensor]*(TS-sat[Wsensor])

				for ip in range(ipx):
					if(nsd[ip]==1):
						if( (sca+pre)==0 and WW[ip]>remainW ):
							WW[ip] = remainW	# remained water
						BW[ip] = WW[ip]
						BS[ip] = SS[ip]
						bdph[ip] = Scover[ip]
						if( pre>0 and T[ip]>0 ):
							TS = diffM*T[ip] + (1-diffM)*BT[ip]
							Qr_plus = BF[ip]*200*TS \
										+ (1-BF[ip])*htr[ip]*(TS-sat[ip])
							Qr = Qr + (Qr_plus)*AS[ip]

					BT[ip] = T[ip]
	
		if( Lday!=Dend or Lyear!=lyear ):	continue

		sim().net.close()
		sim().op.close()
		sim().run.close()
		sim().snow.close()
		sim().data.close()

		RT = RT * interval
		snow_minusT = snow_minusT * interval
		wet_minusT  = wet_minusT * interval
		Qsup = Qsup * interval
		Qr = Qr * interval
		Qr2 = Qr2 * interval

		for i in range(height+1):
			Theight.append( nsca[i] * interval )
		for i in range(slevel+1):
			Tcover.append( nsca[i] * interval )
	
		Tsnow_on  = ntime[0][0] * interval
		Tsnow_off = ntime[0][1] * interval
		Twet_on   = ntime[1][0] * interval
		Twet_off  = ntime[1][1] * interval
		Tdry_on   = ntime[2][0] * interval
		Tdry_off  = ntime[2][1] * interval

		sim().hist.write(str(RT) + ', ')
		sim().hist.write(str(Qsup/area*4.186) + ', ')
		sim().hist.write(str(Qr/area*4.186) + ', ')
		sim().hist.write(str(Qr2*4.186) + ', ')
		sim().hist.write(str(melt) + ', ')
		sim().hist.write(str(area*1000) + '\n')		#[m^2]
		sim().hist.write(str(Tsnow_on) + ', ')
		sim().hist.write(str(Tsnow_off) + ', ')
		sim().hist.write(str(Twet_on) + ', ')
		sim().hist.write(str(Twet_off) + ', ')
		sim().hist.write(str(Tdry_on) + ', ')
		sim().hist.write(str(Tdry_off) + ', ')
		sim().hist.write(str(snow_minusT) + ', ')
		sim().hist.write(str(wet_minusT) + '\n')
		for i in range(height+1):
			sim().hist.write(str(Theight[i]) + ', ')
		sim().hist.write('\n')
		for i in range(slevel+1):
			sim().hist.write(str(Tcover[i]) + ', ')
		sim().hist.write('\n' + str(interval) + ', ')
		sim().hist.write(str(outputStep) + ', ')
		sim().hist.write(str(Dstart) + ', ')
		sim().hist.write(str(lyear) + ', ')
		sim().hist.write(str(Dend) + ', ')
		sim().hist.write('3\n')
		sim().hist.write(str(remainW) + ', ')
		sim().hist.write(str(maxPene) + ', ')
		sim().hist.write(str(abrate0) + '\n')
		sim().hist.write(str(Qs/area/0.86) + ', ')
		sim().hist.write(str(maxT) + ', ')
		sim().hist.write(str(Tsensor) + ', ')
		sim().hist.write(str(Wsensor) + ', ')
		sim().hist.write(str(timer) + '\n')
		sim().hist.write(str(snowT) + ', ')
		sim().hist.write(str(B1) + ', ')
		sim().hist.write(str(stopT1) + ', ')
		sim().hist.write(str(wetT) + ', ')
		sim().hist.write(str(B2) + ', ')
		sim().hist.write(str(stopT2) + ', ')
		sim().hist.write(str(dryT) + ', ')
		sim().hist.write(str(B3) + ', ')
		sim().hist.write(str(stopT3) + ', ')
		sim().hist.write(str(offT) + '\n')

		for ip in range(ipx):
			sim().net2.write(str(ip) + ', ')
			sim().net2.write(str(BT[ip]) + ', ')
			sim().net2.write(str(ID[ip]) + ', ')
			sim().net2.write(str(C[ip]) + ', ')
			sim().net2.write(str(CL[ip]) + '\n')
			sim().net2.write(str(nsd[ip]) + ', ')
			sim().net2.write(str(AS[ip]) + ', ')
			sim().net2.write(str(nqd[ip]) + ', ')
			sim().net2.write(str(QQ[ip]) + '\n')
			sim().net2.write(str(iopt[ip]) + ', ')
			sim().net2.write(str(iopq[ip]) + ', ')
			sim().net2.write(str(iope[ip]) + '\n')
			sim().net2.write(str(npn[ip]) + '\n')
			for j in range(npn[ip]):
				sim().net2.write(str(NP[ip][j]) + ', ')
				sim().net2.write(str(nkr[ip][j]) + ', ')
				sim().net2.write(str(HR[ip][j]) + ', ')
				sim().net2.write(str(ndmy) + '\n') 
		sim().net2.write('lpmx = ' + str(lpmx) + '\n')
