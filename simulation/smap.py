#!/usr/bin/env python
# coding : utf-8

from __future__ import print_function
import math
import numpy as np

RD = 3.1416/180
pido = 43.04
pked = 141.33
pcth = 0
pctl = 0.2
tvm = (pcth + pctl)/2
cth = 1.67


class sim:

	def __init__(self):
		self.t = []
		self.bt = []
		self.Q = []
		self.E = []
		self.tset = []
		self.ID = []
		self.QQ = []
		self.npn = []
		self.nsd = []
		self.nqd = []
		self.C = []
		self.CR = []
		self.AS = []
		self.iopt = []
		self.iope = []
		self.iopq = []
		self.CL = []
		self.ict = []
		self.TD = []
		self.QD = []
		self.ED = []
		self.lms = []
		self.lds = []
		self.jdn = []
		self.jtd = []
		self.jqd = []
		self.jed = []
		self.QE = []
		self.BF = []
		self.htr = []
		self.NP = []
		self.nkr = []
		self.HR = []
		self.S1 = []
		self.S2 = []
		self.TT = []
		self.ER = []
		self.flr = []
		self.IL = []
		self.BW = []
		self.BS = []
		self.SS = []
		self.WW = []
		self.EV = []
		self.dph = []
		self.bdph = []
		self.mlt = []
		self.sat = []
		self.bst = []
		self.NC = []
		self.smf = []
		self.itmp = []
		self.islr = []
		self.ipre = []
		self.IX0 = []
		self.irdv = []
		self.irfv = []
		self.iunr = []
		self.ifko = []
		self.ivel = []
		self.cfd = []
		self.csca = []
		self.nfd = []
		self.nsca = []
		self.ntime = []
		self.dpd = []
		self.ssd = []
		self.tsd = []

		fn1 = 'file.net'
		self.net = open(fn1, 'r')
		fn2 = 'op.csv'
		self.op = open(fn2, 'r')
		fn3 = 'smap.run'
		self.run = open(fn3, 'r')
		fn4 = 'sekisetu.csv'
		self.sekisetu = open(fn4, 'r')


if __name__ == '__main__':

	sse = 0
	RT = 0
	ifdx = 0
	iscax = 0
	etsx = 0
	QF = 0
	qf2 = 0
	idr = 2
	irt = 0
	ssmft = 0

	cfd = np.zeros(100)
	csca = np.zeros(100)

	ctime = np.zeros((2, 5))


	smap1 = sim().run.readline().split('\t')
	dtm = float(smap1[0])
	lpr0 = float(smap1[1])
	CF = float(smap1[2])
	CK = float(smap1[3])
	lpsx = float(smap1[4])
	print(dtm, lpr0, CF, CK, lpsx)

	idx = int( 1/(dtm+0.4) )
	lpr = 1.0
	dlt = 0.0001

	smap2 = sim().run.readlines()[1].split('\t')
	lst = float(smap2[0])
	lye = float(smap2[1])
	lend = float(smap2[2])
	print(lst, lye, lend)

	ND = int( sim().run.readlines()[2] )
	print(ND)
	for i in range(ND):
		smaps = sim().run.readlines()[3+i].split('\t')
	
	LD = lst - 1
	lyr = 0

	smap3 = sim().run.readlines()[3+ND].split('\t')
	mode = int(smap3[0])
	fname = smap3[1]

	smap4 = sim().run.readlines()[4+ND].split('\t')
	QB = smap4[0]
	WF = smap4[1]
	temp = smap4[2]
	I1 = smap4[3]
	I2 = smap4[4]
	tmoff = smap4[5]
	iroff = smap4[6]
	irset = smap4[7]

	smap5 = sim().run.readlines()[5+ND].split('\t')
	A1 = smap4[0]
	B1 = smap4[1]
	C1 = smap4[2]
	A2 = smap4[3]
	B2 = smap4[4]
	C2 = smap4[5]
	A3 = smap4[6]
	B3 = smap4[7]
	C3 = smap4[8]
	taoff = smap4[9]

	smap6 = sim().run.readlines()[6+ND].split('\t')
	flo = smap4[0]
	hix = smap4[1]
	ar0 = smap4[2]
	cvel = smap4[3]
	dph0 = smap4[4]
	ymid = smap4[5]

	if(mode==1):
		sim().sekisetu
